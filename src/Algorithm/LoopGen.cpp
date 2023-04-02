#include "LoopGen.h"
#include <omp.h>

#define COMPUTE_NEW_PLANELOOP 1
#define COMPUTE_NEW_ENERGY 1
#define PRINT_WHY_EXIT 0
namespace LoopGen
{
#pragma region initialization
	void LoopGen::GetPositionFromLoop(const std::vector<VertexLayer*>& loop, Eigen::VectorXd xyz[3])
	{
		int n = loop.size() - 1;
		xyz[0].resize(n); xyz[1].resize(n); xyz[2].resize(n);
		for (int i = 0; i < n; ++i)
		{
			auto& pos = mesh->point(loop[i]->v);
			xyz[0](i) = pos[0];
			xyz[1](i) = pos[1];
			xyz[2](i) = pos[2];
		}
	}

	double LoopGen::EvaluateLoopSimilarity(Eigen::Matrix3Xd& loop0, Eigen::Matrix3Xd& loop1, double u, int begin_seg)
	{
		//对两个loop重新采样，对比采样点的切向，从而定义相似性
		int n = loop0.cols() + loop1.cols();
		auto loopLength = [&](Eigen::Matrix3Xd& loop, Eigen::VectorXd& seg, int mark, Eigen::Vector3d &barycenter)
		{
			int cols = loop.cols();
			double sum = 0;
			seg.resize(cols);
			barycenter.setZero();
			for (int i = 0; i < cols; ++i)
			{
				seg(i) = (loop.col((mark + i + 1) % cols) - loop.col((mark + i) % cols)).norm();
				sum += seg(i);
				barycenter += loop.col(i);
			}
			barycenter /= cols;
			return sum;
		};
		auto assembleSampling = [&](Eigen::Matrix3Xd &sampling, double u0, int mark, Eigen::Matrix3Xd &loop)
		{
			int cols = loop.cols();
			Eigen::VectorXd seg;
			Eigen::Vector3d barycenter; 
			double len = loopLength(loop, seg, mark, barycenter);
			double step = len / n;
			Eigen::Vector3d vec = (loop.col((mark + 1) % cols) - loop.col(mark % cols)).normalized();
			int r = 0;
			double l = seg(0);
			sampling.resize(3, n);
			sampling.col(0) = (loop.col(mark % cols) + u0 * vec - barycenter) / len;
			for (int i = 1; i < n; ++i)
			{
				u0 += step;
				if (u0 > l)
				{
					while (u0 > l)
					{
						++r; ++mark;
						l += seg(r % cols);
					}
					vec = (loop.col((mark + 1) % cols) - loop.col(mark % cols)).normalized();
				}
				//sampling.col(i) = (loop.col(mark % cols) + u0 * vec - source).normalized();
				//sampling.col(i) = (loop.col(mark % cols) + u0 * vec - barycenter) / len;
				sampling.col(i) = (loop.col((mark + 1) % cols) + (u0 - l)*vec - barycenter) / len;
			}
		};
		Eigen::Matrix3Xd sampling0, sampling1;
		assembleSampling(sampling0, 0, 0, loop0);
		assembleSampling(sampling1, u, begin_seg, loop1);
		return ICP_Energy(sampling0, sampling1);
	}

	void LoopGen::LeastSquarePlane(Eigen::VectorXd xyz[3], double plane[4])
	{
		//xyz中第一个点是loop的出发点，这里加强其权重使平面能与之接近，最后再移动平面使之穿过出发点
		int n = xyz[0].size();
		//position.col(0) *= n;
		auto deviation = [&](int k)
		{
			return (xyz[k].array() - xyz[k].sum() / n).abs().sum();
		};
		double dev[3] = { deviation(0),deviation(1),deviation(2) };
		int minId = dev[0] < dev[1] ? (dev[0] < dev[2] ? 0 : 2) : (dev[1] < dev[2] ? 1 : 2);
		auto& p0 = xyz[minId];
		auto& p1 = xyz[(minId + 1) % 3];
		auto& p2 = xyz[(minId + 2) % 3];
		Eigen::Matrix3d m;
		m.col(0) << p1.dot(p1), p1.dot(p2), p1.sum();
		m.col(1) << m(1, 0), p2.dot(p2), p2.sum();
		m.col(2) << m(2, 0), m(2, 1), n;
		Eigen::Vector3d right(p0.dot(p1), p0.dot(p2), p0.sum());
		//加强权重
		if (1)
		{
			n *= n;
			m(0, 0) += n * p1(0) * p1(0); m(0, 1) += n * p1(0) * p2(0); m(0, 2) += n * p1(0);
			m(1, 0) = m(0, 1); m(1, 1) += n * p2(0) * p2(0); m(1, 2) += n * p2(0);
			m(2, 0) = m(0, 2); m(2, 1) = m(1, 2); m(2, 2) += n;
			right(0) += n * p1(0) * p0(0); right(1) += n * p2(0) * p0(0); right(2) += n * p0(0);
		}
		right = m.ldlt().solve(-right);
		double invnorm = 1.0 / sqrt(right(0) * right(0) + right(1) * right(1) + 1);
		right *= invnorm;
		switch (minId)
		{
		case 0:
			plane[0] = invnorm; plane[1] = right(0); plane[2] = right(1);
			break;
		case 1:
			plane[0] = right(1); plane[1] = invnorm; plane[2] = right(0);
			break;
		case 2:
			plane[0] = right(0); plane[1] = right(1); plane[2] = invnorm;
			break;
		}
		//plane[3] = right(2);
		//移动平面
		plane[3] = -(plane[0] * xyz[0](0) + plane[1] * xyz[1](0) + plane[2] * xyz[2](0));
	}

	double LoopGen::EvaluatePlanarity(Eigen::VectorXd xyz[3], double plane[4])
	{
		int n = xyz[0].size();
		double sum = 0;
		for (int i = 0; i < n; ++i)
		{
			sum += fabs(xyz[0](i) * plane[0] + xyz[1](i) * plane[1] + xyz[2](i) * plane[2] + plane[3]);
		}
		return sum / n;
	}

	double LoopGen::RefineLoopByPlanarity(std::vector<VertexLayer*>& loop, PlaneLoop& planar_loop)
	{
		planar_loop.clear();
		Eigen::VectorXd xyz[3];
		GetPositionFromLoop(loop, xyz);
		double plane[4];
		LeastSquarePlane(xyz, plane);

		auto dis = [&](Vec3d& pos)
		{
			return plane[0] * pos[0] + plane[1] * pos[1] + plane[2] * pos[2] + plane[3];
		};

		//从起始点出发，不断加入端点落在平面两侧的边
		HalfedgeLayer* hl_begin = loop[0]->hl;
		HalfedgeLayer* hl_transfer = hl_begin;
		double s0 = dis(mesh->point(m4.verticelayers[hl_transfer->next->to].v));
		double distance[2];
		PointOnHalfedgeLayer pohl[2];
		int id[2] = { 0,0 };
		do
		{
			double s1 = dis(mesh->point(m4.verticelayers[hl_transfer->to].v));
			if (s0*s1 < 0)
			{
				if (s0 > 0)
				{
					pohl[0].hl = hl_transfer->next;
					pohl[0].c = s0 / (s0 - s1);
					distance[0] = s0; distance[1] = s1;
					++id[0];
				}
				else
				{
					pohl[1].hl = hl_transfer->next->oppo;
					pohl[1].c = s1 / (s1 - s0);
					++id[1];
				}
			}
			hl_transfer = hl_transfer->oppo->next;
			s0 = s1;
		} while (hl_transfer != hl_begin);

		if (id[0] != 1 || id[1] != 1)
		{
			return YYSS_INFINITE;
		}

		//检查出发点到pohl[0]的方向与搜索loop的方向是否相同，若不相同，则调换pohl[0]和pohl[1]
		auto &v0 = cf->getCrossField().col(pohl[0].hl->left);
		auto v1 = pohl[0].point(m4) - mesh->point(loop.front()->v);
		if (v0(0)*v1[0] + v0(1)*v1[1] + v0(2)*v1[2] < 0)
		{
			std::swap(pohl[0], pohl[1]);
			for (int i = 0; i < 4; ++i) 
				plane[i] *= -1.0;
			for (int i = 0; i < 2; ++i)
			{
				pohl[i].hl = pohl[i].hl->oppo;
				pohl[i].c = 1 - pohl[i].c;
			}
			distance[0] = dis(mesh->point(m4.verticelayers[pohl[0].hl->to].v));
			distance[1] = dis(mesh->point(m4.verticelayers[pohl[0].hl->from].v));
		}

		planar_loop.emplace_back(pohl[0].hl, pohl[0].c);
		hl_transfer = pohl[0].hl;
		while (hl_transfer->h != pohl[1].hl->h)
		{
			double s = dis(mesh->point(m4.verticelayers[hl_transfer->oppo->prev->from].v));
			if (s > 0)
			{
				hl_transfer = hl_transfer->oppo->next;
				distance[0] = s;
			}
			else
			{
				hl_transfer = hl_transfer->oppo->prev;
				distance[1] = s;
			}
			planar_loop.emplace_back(hl_transfer, distance[0] / (distance[0] - distance[1]));
		}
		if (hl_transfer != pohl[1].hl)
		{
			planar_loop.clear();
			return YYSS_INFINITE;
		}
		return EvaluatePlanarity(xyz, plane);
	}

	bool LoopGen::FieldAligned_PlanarLoop(/*M4 &m4, */VertexLayer* vl, std::vector<VertexLayer*> &path, BoolVector &break_info)
	{
		layernode_pq pq;

		int nvl = m4.verticelayers.size();
		int vid;
		std::vector<double> distance(nvl, YYSS_INFINITE);
		std::vector<int> count(nvl, 0);
		std::vector<HalfedgeLayer*> prev(nvl, nullptr);
		BoolVector visited(nvl, false);

		HalfedgeLayer* hl_begin = vl->hl;
		HalfedgeLayer* hl_transfer = hl_begin;
		//Eigen::Vector4d plane_func; plane_func.setZero();
		auto &crossfield = cf->getCrossField();
		do
		{
			double w = m4.weight(hl_transfer->id);
			if (w < YYSS_INFINITE)
			{
				vid = m4.verticelayers[hl_transfer->to].v.idx();
				if (!m4.sing_flag[vid])
				{
					int toid = hl_transfer->to;
					pq.emplace(toid, w, ++count[toid]);
					distance[toid] = w;
					prev[toid] = hl_transfer;
					//vid = m4.verticemap[vid];
					vid *= 4;
					for (int i = 0; i < 4; ++i)
					{
						if (vid + i == toid)
							continue;
						visited[vid + i] = true;
					}
				}
			}
			//plane_func.head(3) += crossfield.col(m4.conj_hl(hl_transfer, 1)->left);
			hl_transfer = hl_transfer->prev->oppo;
		} while (hl_transfer != hl_begin);
		//plane_func.head(3).normalize();
		OpenMesh::Vec3d pos = mesh->point(vl->v);
		//plane_func(3) = -(plane_func(0)*pos[0] + plane_func(1)*pos[1] + plane_func(2)*pos[2]);
		/*auto plane_dist = [&](OpenMesh::Vec3d &p)
		{
			return fabs(plane_func(0)*p[0] + plane_func(1)*p[1] + plane_func(2)*p[2] + plane_func(3));
		};*/
		VertexLayer* last_vl = nullptr;
		while (true)
		{
			LayerNode ln;
			do
			{
				if (pq.empty())
					return false;
				ln = pq.top();
				pq.pop();
			} while (ln.count != count[ln.id]);
			int fromid = ln.id;
			//if (fromid == vl->id)
			if (break_info[fromid])
			{
				last_vl = &m4.verticelayers[fromid];
				break;
			}

			hl_begin = m4.verticelayers[fromid].hl;
			hl_transfer = hl_begin;
			do
			{
				double w = m4.weight(hl_transfer->id);
				if (w < YYSS_INFINITE)
				{
					vid = m4.verticelayers[hl_transfer->to].v.idx();
					if (!m4.sing_flag[vid] && !visited[hl_transfer->to])
					{
						int toid = hl_transfer->to;
						/*double ratio = plane_dist(mesh->point(mesh->vertex_handle(vid))) / avg_len;
						if (ratio > 2.0)
						{
							w *= ratio / 2;
						}*/
						if (distance[fromid] + w < distance[toid])
						{
							distance[toid] = distance[fromid] + w;
							pq.emplace(toid, distance[toid], ++count[toid]);
							prev[toid] = hl_transfer;
							//vid = m4.verticemap[vid];
							vid *= 4;
							for (int i = 0; i < 4; ++i)
							{
								if (vid + i == toid)
									continue;
								visited[vid + i] = true;
							}
						}
					}
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}

		path.clear();
		path.push_back(last_vl);
		VertexLayer* prev_layer = &m4.verticelayers[prev[last_vl->id]->from];
		while (prev_layer != vl)
		{
			path.push_back(prev_layer);
			prev_layer = &m4.verticelayers[prev[prev_layer->id]->from];
		}
		path.push_back(vl);
		std::reverse(path.begin(), path.end());
		return true;
	}

	void LoopGen::InitializeField()
	{
#if 1
		std::ifstream file_reader;
		std::string field_file = "../resource//field//" + model_name + ".field";
		file_reader.open(field_file, std::ios::in);
		if (file_reader.good())
		{
			file_reader.close();
			cf = new crossField(mesh, field_file);
			cf->read_field();
			cf->initFieldInfo();
		}
		else
		{
			file_reader.close();
			cf = new crossField(mesh, field_file);
			cf->initMeshInfo();
			cf->setCurvatureConstraint();
			cf->setField();
			cf->initFieldInfo();
			cf->write_field();
		}
#else
		cf = new crossField(mesh, field_file);
		cf->initMeshInfo();
		cf->setCurvatureConstraint();
		cf->setField();
		cf->initFieldInfo();
#endif
		dprint("Initialize Field Done!");
	}

	void LoopGen::InitializePlaneLoop()
	{
		m4.set_base(mesh, cf);
		m4.init();
		m4.update();
		m4.set_weight();
		InfoOnMesh.resize(m4.verticelayers.size());
		int count = 0;
		for (int i = 0; i < InfoOnMesh.size(); ++i)
		{
			InfoOnMesh[i].id = i;
			InfoOnMesh[i].energy = YYSS_INFINITE;
			switch (i % 4)
			{
			case 0:
				InfoOnMesh[i].dir = Forward;
				InfoOnMesh[i].plid = count;
				++count;
				break;
			case 1:
				InfoOnMesh[i].dir = Forward;
				InfoOnMesh[i].plid = count;
				break;
			case 2:
				InfoOnMesh[i].dir = Reverse;
				InfoOnMesh[i].plid = count - 1;
				break;
			case 3:
				InfoOnMesh[i].dir = Reverse;
				InfoOnMesh[i].plid = count;
				++count;
				break;
			}
		}
		int nv = mesh->n_vertices();

		//pls.resize((m4.verticelayers.size() - cf->getSingularity().size()) / 2 + cf->getSingularity().size());
		pls.resize(m4.verticelayers.size() / 2);

		tr.tog();
		if (COMPUTE_NEW_PLANELOOP || !ReadPlaneLoop(m4, pls, model_name, mesh))
		{
#pragma omp parallel for
			for (int i = 0; i < nv; ++i)
			{
				if (m4.sing_flag[i])
					continue;
				//VertexLayer* vl = &m4.verticelayers[m4.verticemap[i]];
				VertexLayer* vl = &m4.verticelayers[i * 4];
				std::vector<VertexLayer*> loop;
				BoolVector break_info(m4.verticelayers.size(), false);
				break_info[vl->id] = true;
				if (FieldAligned_PlanarLoop(vl, loop, break_info))
					RefineLoopByPlanarity(loop, pls[InfoOnMesh[vl->id].plid]);
				break_info[vl->id] = false;
				++vl;
				break_info[vl->id] = true;
				if (FieldAligned_PlanarLoop(vl, loop, break_info))
					RefineLoopByPlanarity(loop, pls[InfoOnMesh[vl->id].plid]);
			}
			WritePlaneLoop(pls, model_name, mesh);
		}
		tr.out("Time of Initializing Planar Loops on All Vertices:");
	}

	void LoopGen::InitializeSimilarityEnergy(bool re_compute)
	{
		auto assembleLoop = [&](Vec3d &start, PlaneLoop &pl, DIRECTION dir, Eigen::Matrix3Xd & loop)
		{
			loop.resize(3, 1 + pl.size());
			loop.col(0) << start[0], start[1], start[2];
			int c = 0;
			if (dir == Forward)
			{
				for (auto itr = pl.begin(); itr != pl.end(); ++itr)
				{
					auto pos = itr->point(m4);
					loop.col(++c) << pos[0], pos[1], pos[2];
				}
			}
			else
			{
				for (auto itr = pl.rbegin(); itr != pl.rend(); ++itr)
				{
					auto pos = itr->point(m4);
					loop.col(++c) << pos[0], pos[1], pos[2];
				}
			}
		};

		tr.tog();
		similarity_energy.resize(m4.halfedgelayers.size(), YYSS_INFINITE);
		if (re_compute || COMPUTE_NEW_ENERGY || !ReadEnergy(similarity_energy, model_name))
		{
#pragma omp parallel for
			for (int k = 0; k < mesh->n_edges(); ++k)
			{
				for (int m = 0; m < 2; ++m)
				{
					int hlid = k * 8 + m;
					auto& hl = m4.halfedgelayers[hlid];
					auto& fiov = InfoOnMesh[hl.from];
					auto& tiov = InfoOnMesh[hl.to];
					if (pls[fiov.plid].empty() || pls[tiov.plid].empty())
					{
						continue;
					}
					Eigen::Matrix3Xd loop0, loop1;
					assembleLoop(mesh->point(m4.verticelayers[hl.from].v), pls[fiov.plid], fiov.dir, loop0);
					assembleLoop(mesh->point(m4.verticelayers[hl.to].v), pls[tiov.plid], tiov.dir, loop1);
					auto& proj_pos = loop0.col(0);
					int id = -1;
					int cols = loop1.cols();
					double u0 = 0;
					for (int j = 0; j < 2; ++j)
					{
						Eigen::Vector3d pos0 = loop1.col(j + 1) - loop1.col(j);
						double dot0 = pos0.dot(proj_pos - loop1.col(j));
						double dot1 = pos0.dot(proj_pos - loop1.col(j + 1));
						if (dot0 > 0 && dot1 < 0)
						{
							id = j;
							u0 = dot0 / (dot0 - dot1) * pos0.norm();
							break;
						}
						pos0 = loop1.col((cols - j) % cols) - loop1.col(cols - j - 1);
						dot0 = pos0.dot(proj_pos - loop1.col(cols - j - 1));
						dot1 = pos0.dot(proj_pos - loop1.col((cols - j) % cols));
						if (dot0 > 0 && dot1 < 0)
						{
							id = cols - j - 1;
							u0 = dot0 / (dot0 - dot1) * pos0.norm();
							break;
						}
					}
					if (id == -1)
						id = 0;
					similarity_energy[hlid] = EvaluateLoopSimilarity(loop0, loop1, u0, id);
					similarity_energy[hlid + 2] = similarity_energy[hlid];
					similarity_energy[hl.oppo->id] = similarity_energy[hlid];
					similarity_energy[(hl.oppo->id / 4) * 4 + ((hl.oppo->id + 2) % 4)] = similarity_energy[hlid];
				}
			}
			WriteEnergy(similarity_energy, model_name);
		}
		tr.out("time of computing similarity energy:");

		tr.tog();
#if 1
		for (auto &iov : InfoOnMesh)
			iov.energy = 0;
		for (auto &hl : m4.halfedgelayers)
		{
			InfoOnMesh[hl.from].energy += similarity_energy[hl.id];
			InfoOnMesh[hl.to].energy += similarity_energy[hl.id];
		}
		for (auto &iov : InfoOnMesh)
			iov.energy /= 8 * mesh->valence(m4.verticelayers[iov.id].v);
#else
		for (auto &iov : InfoOnMesh)
			iov.energy = YYSS_INFINITE;
		for (auto &hl : m4.halfedgelayers)
		{
			InfoOnMesh[hl.from].energy = std::min(InfoOnMesh[hl.from].energy, similarity_energy[hl.id]);
			InfoOnMesh[hl.to].energy = std::min(InfoOnMesh[hl.to].energy, similarity_energy[hl.id]);
		}
#endif
		tr.out("time of setting vertex energy:");
	}
#pragma endregion 

	void LoopGen::ConstructInitialRegion(InfoOnVertex* iov, LocalParametrization &lp)
	{
		int nvl = m4.verticelayers.size();
		int nfl = m4.facelayers.size();
		auto& pl = pls[iov->plid];
		std::vector<std::vector<int>> advancing_front[2];
		std::vector<int> hierarchy_vertex[2];
		hierarchy_vertex[0].reserve(pl.size() + 3);
		hierarchy_vertex[1].reserve(pl.size() + 3);
		BoolVector visited_vl(m4.verticelayers.size(), false);
		visited_vl[iov->id] = true;
		for (auto& pohl : pl)
		{
			if (!visited_vl[pohl.hl->from])
			{
				hierarchy_vertex[0].push_back(pohl.hl->from);
				visited_vl[pohl.hl->from] = true;
			}
			if (!visited_vl[pohl.hl->to])
			{
				hierarchy_vertex[1].push_back(pohl.hl->to);
				visited_vl[pohl.hl->to] = true;
			}
		}
		advancing_front[0].push_back(std::move(hierarchy_vertex[0]));
		advancing_front[1].push_back(std::move(hierarchy_vertex[1]));
		while (true)
		{
			const auto& current_af = advancing_front[0].back();
			std::vector<int> hierarchy;
			bool exceed_energy = false;
			for (auto caf : current_af)
			{
				VertexLayer* vl = &m4.verticelayers[caf];
				HalfedgeLayer* hl_begin = vl->hl;
				HalfedgeLayer* hl_transfer = hl_begin;
				do
				{
					int toid = hl_transfer->to;
					if (!visited_vl[toid])
					{
						if (m4.sing_flag[toid / 4])
							goto target0;
						if (InfoOnMesh[toid].energy > energy_threshold)
							exceed_energy = true;
						visited_vl[toid] = true;
						hierarchy.push_back(toid);
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);
			}
			if (hierarchy.empty())
				break;
			advancing_front[0].push_back(std::move(hierarchy));
			if (exceed_energy)
				break;
		}
	target0:;
		while (true)
		{
			const auto& current_af = advancing_front[1].back();
			std::vector<int> hierarchy;
			bool exceed_energy = false;
			for (auto caf : current_af)
			{
				VertexLayer* vl = &m4.verticelayers[caf];
				HalfedgeLayer* hl_begin = vl->hl;
				HalfedgeLayer* hl_transfer = hl_begin;
				do
				{
					int toid = hl_transfer->to;
					if (!visited_vl[toid])
					{
						if (m4.sing_flag[toid / 4])
							goto target1;
						if (InfoOnMesh[toid].energy > energy_threshold)
							exceed_energy = true;
						visited_vl[toid] = true;
						hierarchy.push_back(toid);
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);
			}
			if (hierarchy.empty())
				break;
			advancing_front[1].push_back(std::move(hierarchy));
			if (exceed_energy)
				break;
		}
	target1:;

		auto& new_vertex = lp.new_vertex;
		auto& new_face = lp.new_face;
		auto& newv_flag = lp.new_v_flag; newv_flag.resize(nvl, false);
		auto& newf_flag = lp.new_f_flag; newf_flag.resize(nfl, false);
		auto& regionv_flag = lp.region_v_flag;
		auto& regionf_flag = lp.region_f_flag;
		for (int i = 0; i < 2; ++i)
		{
			for (auto& ss : advancing_front[i])
			{
				for (auto& tt : ss)
				{
					new_vertex.push_back(&m4.verticelayers[tt]);
					newv_flag[tt] = true;
				}
			}
		}

		for (auto new_v : new_vertex)
		{
			HalfedgeLayer* hl_begin = new_v->hl;
			HalfedgeLayer* hl_transfer = hl_begin;
			do
			{
				int vf_id = hl_transfer->left;
				if (!(regionf_flag[vf_id] || newf_flag[vf_id]))
				{
					HalfedgeLayer* hb = m4.facelayers[vf_id].hl;
					HalfedgeLayer* ht = hb;
					do
					{
						if (!(regionv_flag[ht->to] || newv_flag[ht->to]))
							goto target2;
						ht = ht->next;
					} while (ht != hb);
					new_face.push_back(&m4.facelayers[vf_id]);
					newf_flag[vf_id] = true;
				}
			target2:;
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}


		auto& x_axis = lp.x_axis;
		auto& y_axis = lp.y_axis;
		auto& crossfield = cf->getCrossField();
		for (auto fl : new_face)
		{
			int flid = fl->id;
			int fid = m4.facelayers[flid].f.idx();
			x_axis.col(fid) = crossfield.col(flid);
			y_axis.col(fid) = crossfield.col(fid * 4 + (flid + 1) % 4);
		}

		//看情况决定是否优化场
		if (cset.has_face.size() == 0)
			return;
		std::vector<FaceLayer*> opt_face;
		BoolVector opt_flag(m4.facelayers.size(), false);
		BoolVector constraint_flag(m4.facelayers.size(), false);
		for (auto fl : new_face)
		{
			if (cset.has_face[fl->f.idx()])
			{
				constraint_flag[fl->id] = true;
			}
			else
			{
				opt_face.push_back(fl);
				opt_flag[fl->id] = true;
			}
		}
		if (!opt_face.empty())
			ResetLocalField(lp, opt_face, opt_flag, constraint_flag);
	}

	void LoopGen::AssembleSampling(VertexLayer* vl, Eigen::Matrix3Xd &sampling, LocalParametrization &lp, int loop_sampling_num)
	{
		double step = 1.0 / loop_sampling_num;
		auto setData = [&](double t0, double t1, double c, Eigen::Matrix3Xd& sampling, double fromu, double& tou, OpenMesh::Vec3d& frompos, OpenMesh::Vec3d& topos)
		{
			if (fabs(t0 - t1) > 0.5) { if (t0 < t1) { t0 += 1.0; } else { t1 += 1.0; } }
			tou = c * t0 + (1 - c) * t1;
			tou -= std::floor(tou);
			if (fromu > tou && fabs(fromu - tou) < 0.5) return;

			int u0 = std::floor(fromu * loop_sampling_num);
			int u1 = std::floor(tou * loop_sampling_num);
			if (u0 == u1) return;
			//OpenMesh::Vec3d ev = (topos - frompos).normalized();
			if (u0 < u1) for (int i = u0 + 1; i <= u1; ++i) {
				if (fabs(fromu - tou) < YYSS_FAIRLY_SMALL)
					c = 0.5;
				else
					c = (i*step - fromu) / (tou - fromu);
				for (int j = 0; j < 3; ++j)
					sampling(j, i) = (1 - c)*frompos[j] + c * topos[j];
			}
			else {
				for (int i = u0 + 1; i < loop_sampling_num; ++i) {
					if (fabs(1.0 - fromu) < YYSS_FAIRLY_SMALL)
						c = 0.5;
					else
						c = (i*step - fromu) / (1.0 - fromu);
					for (int j = 0; j < 3; ++j)
						sampling(j, i) = (1 - c)*frompos[j] + c * topos[j];
				}
				for (int i = 0; i <= u1; ++i) {
					if (u1 == 0)
						c = 0.5;
					else
						c = i * step / tou;
					for (int j = 0; j < 3; ++j)
						sampling(j, i) = (1 - c)*frompos[j] + c * topos[j];
				}
			}
		};
		sampling.resize(3, loop_sampling_num); sampling.setZero();
		double len = 0;
		OpenMesh::Vec3d barycenter(0.0, 0.0, 0.0);
		OpenMesh::Vec3d frompos, topos;
		double fromu, tou;
		frompos = mesh->point(vl->v);
		fromu = lp.GetRegularU(vl->id);
		for (auto& pl : lp.all_pl[vl->id])
		{
			topos = pl.c * mesh->point(m4.verticelayers[pl.hl->from].v) + (1 - pl.c) * mesh->point(m4.verticelayers[pl.hl->to].v);
			barycenter += topos;
			len += (frompos - topos).norm();
			setData(lp.GetU(pl.hl->from), lp.GetU(pl.hl->to), pl.c, sampling, fromu, tou, frompos, topos);
			frompos = topos;
			fromu = tou;
		}
		topos = mesh->point(vl->v);
		barycenter += topos;
		len += (frompos - topos).norm();
		tou = lp.GetRegularU(vl->id);
		setData(tou, tou, 0, sampling, fromu, tou, frompos, topos);

		Eigen::Vector3d bc(barycenter[0], barycenter[1], barycenter[2]);
		bc /= loop_sampling_num;
		for (int i = 0; i < loop_sampling_num; ++i)
			sampling.col(i) = (sampling.col(i) - bc) / len;
	}

	void LoopGen::AssembleSimilarityAngle(VertexLayer* vl, Eigen::VectorXd& sa, LocalParametrization& lp, int loop_fragment_num)
	{
		auto setData = [&](double t0, double t1, double c, Eigen::Matrix3Xd& fragment, double fromu, double& tou, OpenMesh::Vec3d& frompos, OpenMesh::Vec3d& topos)
		{
			if (fabs(t0 - t1) > 0.5) { if (t0 < t1) { t0 += 1.0; } else { t1 += 1.0; } }
			tou = c * t0 + (1 - c) * t1;
			tou -= std::floor(tou);
			if (fromu > tou && fabs(fromu - tou) < 0.5) return;
			int u0 = std::floor(fromu * loop_fragment_num);
			int u1 = std::floor(tou * loop_fragment_num);
			if (u0 == u1) return;
			OpenMesh::Vec3d ev = (topos - frompos).normalized();
			if (u0 < u1) for (int i = u0 + 1; i <= u1; ++i) {
				fragment.col(i) << ev[0], ev[1], ev[2];
			}
			else {
				for (int i = u0 + 1; i < loop_fragment_num; ++i) {
					fragment.col(i) << ev[0], ev[1], ev[2];
				}
				for (int i = 0; i <= u1; ++i) {
					fragment.col(i) << ev[0], ev[1], ev[2];
				}
			}
		};
		Eigen::Matrix3Xd fragment(3, loop_fragment_num); fragment.setZero();
		OpenMesh::Vec3d frompos, topos;
		double fromu, tou;
		frompos = mesh->point(vl->v);
		fromu = lp.GetRegularU(vl->id);
		for (auto& pl : lp.all_pl[vl->id])
		{
			topos = pl.c * mesh->point(m4.verticelayers[pl.hl->from].v) + (1 - pl.c) * mesh->point(m4.verticelayers[pl.hl->to].v);
			setData(lp.GetU(pl.hl->from), lp.GetU(pl.hl->to), pl.c, fragment, fromu, tou, frompos, topos);
			frompos = topos;
			fromu = tou;
		}
		topos = mesh->point(vl->v);
		tou = lp.GetRegularU(vl->id);
		setData(tou, tou, 0, fragment, fromu, tou, frompos, topos);
		sa.resize(loop_fragment_num * (loop_fragment_num + 1) / 2); sa.setZero();
		int count = 0;
		for (int i = 0; i < loop_fragment_num; ++i)
		{
			for (int j = i + 1; j < loop_fragment_num; ++j)
			{
				sa(count++) = conservativeArcCos(fragment.col(i).dot(fragment.col(j)));
			}
		}
	}

	bool LoopGen::SpreadSubRegion(LocalParametrization& lp, bool grow_flag[2])
	{
		int nvl = m4.verticelayers.size();
		int nfl = m4.facelayers.size();

		auto& newv_flag = lp.new_v_flag;
		auto& newf_flag = lp.new_f_flag;
		auto& new_vertex = lp.new_vertex;
		auto& new_face = lp.new_face;
		BoolVector visited_v = lp.region_v_flag;
		BoolVector visited_f = lp.region_f_flag;
		for (auto new_v : new_vertex)
		{
			visited_v[new_v->id] = true;
			newv_flag[new_v->id] = false;
		}
		for (auto new_f : new_face)
		{
			visited_f[new_f->id] = true;
			newf_flag[new_f->id] = false;
		}
		//dprint("new vertex:", lp.new_vertex.size());
		//提取含有完整loop的区域
		if (lp.region_vertex.size() == 1)
			RefineLoopByParametrization(lp.region_vertex.front(), lp, visited_v, visited_f);
		std::vector<VertexLayer*> vertex_cache;
		BoolVector vertex_cache_flag = lp.region_v_flag;

		for (auto new_v : new_vertex)
		{
			int se = lp.all_pl[new_v->id].size();
			if (RefineLoopByParametrization(new_v, lp, visited_v, visited_f))
			{
				vertex_cache_flag[new_v->id] = true;
				vertex_cache.push_back(new_v);
			}
		}

		//dprint("vertex cache:", vertex_cache.size());
		v_cache_flag = vertex_cache_flag;
		v_cache = vertex_cache;
#if PRINT_DEBUG_INFO
		dprint("提取含有完整loop的区域");
#endif
		if (lp.region_face.empty())
		{
			auto hl_begin = lp.region_vertex.front()->hl;
			auto hl_transfer = hl_begin;
			do
			{
				if (!vertex_cache_flag[hl_transfer->to])
				{
#if PRINT_WHY_EXIT
					dprint("因起始点周围缺少点而退出");
#endif
					return false;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			if (!CheckCylinderTopology(vertex_cache, vertex_cache_flag, lp.grow_dir))
			{
#if PRINT_WHY_EXIT
				dprint("因拓扑存在问题而退出");
#endif
				return false;
			}
		}
		if (vertex_cache.empty())
		{
#if PRINT_WHY_EXIT
			dprint("因新增点为空而退出");
#endif
			return false;
		}

#if PRINT_DEBUG_INFO
		dprint("拓扑检查");
#endif

		auto& regionv_flag = lp.region_v_flag;
		auto& regionf_flag = lp.region_f_flag;
		auto& grow_dir = lp.grow_dir;

		//检测相似性能量
		auto& all_pl = lp.all_pl;
#if USE_NEW_SIMILARITY_ENERGY
		auto& normal_sampling = lp.normal_sampling;
		int loop_sampling_num = all_pl[lp.region_vertex.front()->id].size();
		if (vertex_cache_flag[lp.region_vertex.front()->id] && !lp.has_ns)
		{
			lp.has_ns = true;
			AssembleSampling(lp.region_vertex.front(), normal_sampling, lp, loop_sampling_num);
		}
#else
		auto &normal_similarity_angle = lp.normal_similarity_angle;
		int loop_fragment_num = all_pl[lp.region_vertex.front()->id].size();
		if (vertex_cache_flag[lp.region_vertex.front()->id] && !lp.has_nsa)
		{
			lp.has_nsa = true;
			AssembleSimilarityAngle(lp.region_vertex.front(), normal_similarity_angle, lp, loop_fragment_num);
		}
#endif
		

		BoolVector if_similarity_energy_low;
		if (lp.region_vertex.size() > 1)
		{
			if_similarity_energy_low.resize(nvl, false);
			int exceed[2] = { 0,0 };
#if PRINT_DEBUG_INFO
#pragma omp parallel for
#else
			double max_ene = 0;
#endif // !PRINT_DEBUG_INFO
			for (int i = 0; i < vertex_cache.size(); ++i)
			{
				int new_id = vertex_cache[i]->id;
				int grow_id = grow_dir[new_id];
				if (!grow_flag[grow_id])
					continue;
#if USE_NEW_SIMILARITY_ENERGY
				Eigen::Matrix3Xd sampling;
				AssembleSampling(vertex_cache[i], sampling, lp, loop_sampling_num);
				double ey = ICP_Energy(normal_sampling, sampling);
				//dprint(i, ey);
				if (ey < energy_threshold)
					if_similarity_energy_low[new_id] = true;
#if PRINT_DEBUG_INFO
				max_ene = std::max(max_ene, ey);
#endif // !PRINT_DEBUG_INFO
#else
				Eigen::VectorXd similarity_angle;
				AssembleSimilarityAngle(vertex_cache[i], similarity_angle, lp, loop_fragment_num);
				double sum = 0;
				for (int j = 0; j < similarity_angle.size(); ++j) {
					sum += fabs(normal_similarity_angle(j) - similarity_angle(j));
				}
				if (sum < energy_threshold * similarity_angle.size())
				{
					if_similarity_energy_low[new_id] = true;
			    }
#endif
			}
			//dprint("max ene:", max_ene);
#if PRINT_DEBUG_INFO
			dprint("检测相似性能量");
			dprint("最高能量为：", max_ene);
#endif
			for (auto vc : vertex_cache)
			{
				if (!grow_flag[grow_dir[vc->id]])
				{
					if_similarity_energy_low[vc->id] = false;
					continue;
				}
				if (!if_similarity_energy_low[vc->id])
				{
					++exceed[grow_dir[vc->id]];
				}
				if (exceed[grow_dir[vc->id]] > 2)
				{
#if PRINT_WHY_EXIT
					dprint("因能量过高,", grow_dir[vc->id], "号方向停止搜索");
#endif
					grow_flag[grow_dir[vc->id]] = false;
				}
			}
	    }
		else
			if_similarity_energy_low.resize(nvl, true);
#if PRINT_DEBUG_INFO
		dprint(grow_flag[0], grow_flag[1]);
#endif
		/*std::ofstream file_writer;
		file_writer.open("C:\\Users\\123\\Desktop\\vscode compare\\new.txt");
		for (int i = 0; i < nvl; ++i)
			file_writer << i << " " << if_similarity_energy_low[i] << std::endl;
		file_writer.close();*/


		for (int i = 0; i < 2; ++i)
		{
			if (grow_flag[i])
			{
				double grad = LoopLenGrad(vertex_cache, lp, vertex_cache_flag, i);
				//dprint("grad:", i, grad);
				if (grad > PI)
				{
#if PRINT_WHY_EXIT
					dprint("因loop梯度过高，", i, "号方向停止搜索");
#endif
					grow_flag[i] = false;
				}
			}
		}
#if PRINT_DEBUG_INFO
		dprint(grow_flag[0], grow_flag[1]);
#endif

		//更新新区域
		int count = lp.region_vertex.size();
		int begin_ = count;
		auto& u_para = lp.GetU();
		auto& v_para = lp.GetV();
		auto& region_vertex = lp.region_vertex;
		auto& vidmap = lp.vidmap;
		for (auto new_v : vertex_cache)
		{
			if (if_similarity_energy_low[new_v->id])
			{
				region_vertex.push_back(new_v);
				int newid = new_v->id;
				regionv_flag[newid] = true;
				u_para(count) = lp.GetU(newid);
				v_para(count) = lp.GetV(newid);
				vidmap[newid] = count++;
			}
		}
		u_para.conservativeResize(count);
		v_para.conservativeResize(count);

		auto& region_face = lp.region_face;
		for (auto new_f : new_face)
		{
			auto hl_begin = new_f->hl;
			auto hl_transfer = hl_begin;
			do
			{
				if (!regionv_flag[hl_transfer->to])
					goto target1;
				hl_transfer = hl_transfer->next;
			} while (hl_transfer != hl_begin);
			int newid = new_f->id;
			region_face.push_back(new_f);
			regionf_flag[newid] = true;
		target1:;
		}

#if PRINT_DEBUG_INFO
		dprint("更新新区域");
		dprint("当前区域含有的顶点数：", region_vertex.size());
#endif
#if 0
		return false;
#endif

		if (!grow_flag[2])
			grow_flag[0] = false;
		if (!grow_flag[3])
			grow_flag[1] = false;
		//若已超过能量阈值或者遇到接近平面的区域，则退出
		if (!grow_flag[0] && !grow_flag[1])
		{
#if PRINT_WHY_EXIT
			dprint("因两个方向搜索停止而退出");
#endif
			return false;
		}

		//扩展advancing_front
		newv_flag.resize(nvl, false);
		new_vertex.clear();
		for (int i = begin_; i < count; ++i)
		{
			auto rvi = region_vertex[i];
			int growid = grow_dir[rvi->id];
			if (!grow_flag[growid])
				continue;
			auto hl_begin = rvi->hl;
			auto hl_transfer = hl_begin;
			do
			{
				int vvid = hl_transfer->to;
				if (!(regionv_flag[vvid] || newv_flag[vvid]))
				{
					if (m4.sing_flag[m4.verticelayers[vvid].v.idx()])
					{
#if PRINT_WHY_EXIT
						dprint("因遇到奇异点，", growid, "号方向停止搜索");
#endif
						grow_flag[growid + 2] = false;
						goto target2;
					}
					newv_flag[vvid] = true;
					new_vertex.push_back(&m4.verticelayers[vvid]);
					grow_dir[vvid] = growid;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}
	target2:;
		begin_ = 0;
		for (int i = 0; i < extend_layer - 1; ++i)
		{
			int end_ = new_vertex.size();
			for (int j = begin_; j < end_; ++j)
			{
				auto rvj = new_vertex[j];
				int growid = grow_dir[rvj->id];
				auto hl_begin = rvj->hl;
				auto hl_transfer = hl_begin;
				do
				{
					int vvid = hl_transfer->to;
					if (!(regionv_flag[vvid] || newv_flag[vvid]))
					{
						if (m4.sing_flag[m4.verticelayers[vvid].v.idx()])
						{
#if PRINT_WHY_EXIT
							dprint("因遇到奇异点，", growid, "号方向停止搜索");
#endif
							grow_flag[growid + 2] = false;
							goto target3;
						}
						newv_flag[vvid] = true;
						new_vertex.push_back(&m4.verticelayers[vvid]);
						grow_dir[vvid] = growid;
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);
			}
			begin_ = end_;
		}
	target3:;
		newf_flag.resize(nfl, false);
		new_face.clear();
		for (auto new_v : new_vertex)
		{
			int newid = new_v->id;
			auto hl_begin = new_v->hl;
			auto hl_transfer = hl_begin;
			do
			{
				int vfid = hl_transfer->left;
				if (!(regionf_flag[vfid] || newf_flag[vfid]))
				{
					auto hb = m4.facelayers[vfid].hl;
					auto ht = hb;
					do
					{
						if (!regionv_flag[ht->to] && !newv_flag[ht->to])
							goto target4;
						ht = ht->next;
					} while (ht != hb);
					newf_flag[vfid] = true;
					new_face.push_back(&m4.facelayers[vfid]);
				target4:;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}
#if PRINT_DEBUG_INFO
		dprint("扩展advancing_front");
#endif
		
		//优化新区域的场
		ResetLocalField(lp, new_face, newf_flag, regionf_flag);
#if PRINT_DEBUG_INFO
		dprint("优化新区域的场");
#endif
		return true;
	}

	void LoopGen::ResetLocalField(LocalParametrization& lp, std::vector<FaceLayer*>& opt_face, BoolVector& opt_flag, BoolVector& constraint_flag)
	{
		typedef std::complex<double> COMPLEX;
		std::vector<Eigen::Triplet<COMPLEX>> triple;
		int count = 0;
		std::vector<int> fidmap(m4.facelayers.size());
		for (auto fl : opt_face)
		{
			fidmap[fl->id] = count++;
		}
		Eigen::VectorXcd b(3 * opt_face.size()); b.setZero();
		count = 0;
		auto& position = cf->getPosition();
		auto& faceBase = cf->getFaceBase();
		auto& x_axis = lp.x_axis;
		auto& y_axis = lp.y_axis;
		for (auto fl : opt_face)
		{
			int flid = fl->id;
			int fid = fl->f.idx();
			auto hl_begin = fl->hl;
			auto hl_transfer = hl_begin;
			do
			{
				if (!mesh->is_boundary(hl_transfer->oppo->h))
				{
					int gid = mesh->face_handle(hl_transfer->oppo->h).idx();
					int glid = hl_transfer->oppo->left;
					if (constraint_flag[glid] || (opt_flag[glid] && flid > glid))
					{
						auto ev = (position.col(mesh->to_vertex_handle(hl_transfer->h).idx()) -
							position.col(mesh->from_vertex_handle(hl_transfer->h).idx())).normalized();
						COMPLEX e_f = COMPLEX(ev.dot(faceBase.col(fid * 2)), -ev.dot(faceBase.col(fid * 2 + 1)));
						COMPLEX e_g = COMPLEX(ev.dot(faceBase.col(gid * 2)), -ev.dot(faceBase.col(gid * 2 + 1)));
						if (opt_flag[flid])
						{
							triple.emplace_back(count, fidmap[flid], e_f);
						}
						else
						{
							COMPLEX dir = COMPLEX(x_axis.col(fid).dot(faceBase.col(2 * fid)), x_axis.col(fid).dot(faceBase.col(2 * fid + 1)));
							b(count) -= e_f * dir;
						}
						if (opt_flag[glid])
						{
							triple.emplace_back(count, fidmap[glid], -e_g);
						}
						else
						{
							COMPLEX dir = COMPLEX(x_axis.col(gid).dot(faceBase.col(2 * gid)), x_axis.col(gid).dot(faceBase.col(2 * gid + 1)));
							b(count) += e_g * dir;
						}
						++count;
					}
				}
				hl_transfer = hl_transfer->next;
			} while (hl_transfer != hl_begin);
		}

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> solver;
		Eigen::SparseMatrix<COMPLEX> A(count, opt_face.size());
		A.setFromTriplets(triple.begin(), triple.end());
		Eigen::SparseMatrix<COMPLEX> AH = A.adjoint();
		solver.compute(AH*A);
		b = AH * b.head(count);
		b = solver.solve(b);
		for (auto fl : opt_face)
		{
			int fid = fl->f.idx();
			double theta = std::arg(b(fidmap[fl->id]));
			x_axis.col(fid) = faceBase.col(fid * 2) * cos(theta) + faceBase.col(fid * 2 + 1) * sin(theta);
			y_axis.col(fid) = faceBase.col(fid * 2) * cos(theta + 0.5 * PI) + faceBase.col(fid * 2 + 1) * sin(theta + 0.5 * PI);
		}
	}

	bool LoopGen::ConstructRegionCut(VertexLayer* vl, BoolVector& visited, std::vector<VertexLayer*>& cut)
	{
		cut.clear();
		cut.push_back(vl);
		vl = m4.another_layer(vl, 1);
		HalfedgeLayer* hl_begin = vl->hl;
		HalfedgeLayer* hl_transfer = hl_begin;	
		while (true)
		{
			double w = YYSS_INFINITE;
			HalfedgeLayer* hl_mark = nullptr;
			do
			{
				if (m4.weight(hl_transfer->oppo->id) < w)
				{
					w = m4.weight(hl_transfer->oppo->id);
					hl_mark = hl_transfer;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			if (!hl_mark)
				return false;
			if (!visited[m4.another_layer(&m4.verticelayers[hl_mark->to], 3)->id])
				break;
			cut.push_back(m4.another_layer(&m4.verticelayers[hl_mark->to], 3));
			hl_begin = hl_mark->oppo;
			hl_transfer = hl_begin;
		}
		std::reverse(cut.begin(), cut.end());

		hl_begin = vl->hl;
		hl_transfer = hl_begin;
		while (true)
		{
			double w = YYSS_INFINITE;
			HalfedgeLayer* hl_mark = nullptr;
			do
			{
				if (m4.weight(hl_transfer->id) < w)
				{
					w = m4.weight(hl_transfer->id);
					hl_mark = hl_transfer;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			if (!hl_mark)
				return false;
			if (!visited[m4.another_layer(&m4.verticelayers[hl_mark->to], 3)->id])
				break;
			cut.push_back(m4.another_layer(&m4.verticelayers[hl_mark->to], 3));
			hl_begin = hl_mark->oppo;
			hl_transfer = hl_begin;
		}
#if PRINT_DEBUG_INFO
		dprint("计算cut");
#endif
		if (cut.size() < 3)
			return false;
		return true;
	}

	void LoopGen::IterativePQ()
	{
		m4.set_base(mesh, cf);
		m4.init();
		m4.update();
		m4.set_weight();
		InfoOnMesh.resize(m4.verticelayers.size());
		int count = 0;
		for (int i = 0; i < InfoOnMesh.size(); ++i)
		{
			InfoOnMesh[i].id = i;
			InfoOnMesh[i].energy = YYSS_INFINITE;
			switch (i % 4)
			{
			case 0:
				InfoOnMesh[i].dir = Forward;
				InfoOnMesh[i].plid = count;
				++count;
				break;
			case 1:
				InfoOnMesh[i].dir = Forward;
				InfoOnMesh[i].plid = count;
				break;
			case 2:
				InfoOnMesh[i].dir = Reverse;
				InfoOnMesh[i].plid = count - 1;
				break;
			case 3:
				InfoOnMesh[i].dir = Reverse;
				InfoOnMesh[i].plid = count;
				++count;
				break;
			}
		}
		int nv = mesh->n_vertices();

		pls.resize(m4.verticelayers.size() / 2);

//		tr.tog();
//		if (COMPUTE_NEW_PLANELOOP || !ReadPlaneLoop(m4, pls, model_name, mesh))
//		{
//#pragma omp parallel for
//			for (int i = 0; i < nv; ++i)
//			{
//				if (m4.sing_flag[i])
//					continue;
//				//VertexLayer* vl = &m4.verticelayers[m4.verticemap[i]];
//				VertexLayer* vl = &m4.verticelayers[i * 4];
//				std::vector<VertexLayer*> loop;
//				if (FieldAligned_PlanarLoop(vl, loop))
//					RefineLoopByPlanarity(loop, pls[InfoOnMesh[vl->id].plid]);
//				++vl;
//				if (FieldAligned_PlanarLoop(vl, loop))
//					RefineLoopByPlanarity(loop, pls[InfoOnMesh[vl->id].plid]);
//			}
//			WritePlaneLoop(pls, model_name, mesh);
//		}
//		tr.out("Time of Initializing Planar Loops on All Vertices:");

		auto assembleLoop = [&](Vec3d &start, PlaneLoop &pl, DIRECTION dir, Eigen::Matrix3Xd & loop)
		{
			loop.resize(3, 1 + pl.size());
			loop.col(0) << start[0], start[1], start[2];
			int c = 0;
			if (dir == Forward)
			{
				for (auto itr = pl.begin(); itr != pl.end(); ++itr)
				{
					auto pos = itr->point(m4);
					loop.col(++c) << pos[0], pos[1], pos[2];
				}
			}
			else
			{
				for (auto itr = pl.rbegin(); itr != pl.rend(); ++itr)
				{
					auto pos = itr->point(m4);
					loop.col(++c) << pos[0], pos[1], pos[2];
				}
			}
		};

		tr.tog();
		similarity_energy.resize(m4.halfedgelayers.size(), YYSS_INFINITE);
		if (COMPUTE_NEW_ENERGY || !ReadEnergy(similarity_energy, model_name))
		{
#pragma omp parallel for
			for (int k = 0; k < mesh->n_edges(); ++k)
			{
				for (int m = 0; m < 2; ++m)
				{
					int hlid = k * 8 + m;
					auto& hl = m4.halfedgelayers[hlid];
					auto& fiov = InfoOnMesh[hl.from];
					auto& tiov = InfoOnMesh[hl.to];
					if (pls[fiov.plid].empty() || pls[tiov.plid].empty())
					{
						continue;
					}
					Eigen::Matrix3Xd loop0, loop1;
					assembleLoop(mesh->point(m4.verticelayers[hl.from].v), pls[fiov.plid], fiov.dir, loop0);
					assembleLoop(mesh->point(m4.verticelayers[hl.to].v), pls[tiov.plid], tiov.dir, loop1);
					auto& proj_pos = loop0.col(0);
					int id = -1;
					int cols = loop1.cols();
					double u0 = 0;
					for (int j = 0; j < 2; ++j)
					{
						Eigen::Vector3d pos0 = loop1.col(j + 1) - loop1.col(j);
						double dot0 = pos0.dot(proj_pos - loop1.col(j));
						double dot1 = pos0.dot(proj_pos - loop1.col(j + 1));
						if (dot0 > 0 && dot1 < 0)
						{
							id = j;
							u0 = dot0 / (dot0 - dot1) * pos0.norm();
							break;
						}
						pos0 = loop1.col((cols - j) % cols) - loop1.col(cols - j - 1);
						dot0 = pos0.dot(proj_pos - loop1.col(cols - j - 1));
						dot1 = pos0.dot(proj_pos - loop1.col((cols - j) % cols));
						if (dot0 > 0 && dot1 < 0)
						{
							id = cols - j - 1;
							u0 = dot0 / (dot0 - dot1) * pos0.norm();
							break;
						}
					}
					if (id == -1)
						id = 0;
					similarity_energy[hlid] = EvaluateLoopSimilarity(loop0, loop1, u0, id);
					similarity_energy[hlid + 2] = similarity_energy[hlid];
					similarity_energy[hl.oppo->id] = similarity_energy[hlid];
					similarity_energy[(hl.oppo->id / 4) * 4 + ((hl.oppo->id + 2) % 4)] = similarity_energy[hlid];
				}
			}
			WriteEnergy(similarity_energy, model_name);
		}
		tr.out("time of computing similarity energy:");

		tr.tog();
		for (auto &iov : InfoOnMesh)
			iov.energy = 0;
		for (auto &hl : m4.halfedgelayers)
		{
			InfoOnMesh[hl.from].energy += similarity_energy[hl.id];
			InfoOnMesh[hl.to].energy += similarity_energy[hl.id];
		}
		for (auto &iov : InfoOnMesh)
			iov.energy /= 8 * mesh->valence(m4.verticelayers[iov.id].v);
		tr.out("time of setting vertex energy:");
	}

	void LoopGen::ConstructCylinder()
	{
		timeRecorder tr;
		vl_pair_pq pq;
		AssembleIOVLoopEnergy(pq);
		

		std::vector<std::vector<int>> region_index(m4.verticelayers.size());
		BoolVector constraint_flag(mesh->n_faces(), false);
		Eigen::Matrix3Xd constraint_dir(3, mesh->n_faces());
		all_plane_loop.resize(m4.verticelayers.size());

		while (true)
		{
			vl_pair ip;
			ip = pq.top(); pq.pop();
			while (!region_index[ip.vl->id].empty())
			{
				if (pq.empty())
				{
					ip.data = YYSS_INFINITE;
					break;
				}
				ip = pq.top(); pq.pop();
			}
			if (ip.data > energy_threshold)
				break;
			//ip.vl = &m4.verticelayers[43159];
			//if (ip.vl->v.idx() != 4616)
			//	continue;
			if (ip.vl->id == 151453)
			{
				int p = 0;
			}
			LocalParametrization lp(m4, ip.vl);
			ConstructInitialRegion(&InfoOnMesh[ip.vl->id], lp);
			//dprint(ip.vl->id, lp.new_vertex.size());
			bool grow_flag[4] = { true, true, true, true };
			do
			{
				/*static int fe = 0;
				++fe;
				dprint("fe", fe);
				dprint("vertices:", lp.region_vertex.size());*/
				BoolVector visited_v = lp.region_v_flag;
				for (auto& ver : lp.new_vertex)
					visited_v[ver->id] = true;
				if (!ConstructRegionCut(ip.vl, visited_v, lp.cut))
					break;
				lp.run(cf->getNormal());
				if (ip.vl->v.idx() == 6440)
				{
					//break;
				}

			} while (SpreadSubRegion(lp, grow_flag));
			//dprint("vertices:", lp.region_vertex.size());
			if (lp.region_vertex.size() > 1)
				lp.modify_cut();

			old_face_flag = lp.region_f_flag;
			old_vert_flag = lp.region_v_flag;
			new_face_flag = lp.new_f_flag;
			new_vert_flag = lp.new_v_flag;
			uv_para[0] = lp.uv[0];
			uv_para[1] = lp.uv[1];
			vertexidmap = lp.vidmap;
			xaxis = lp.x_axis;
			cut_vertex_flag = lp.cutv_flag;
			growDIR = lp.grow_dir;

			int seed_index = seed_vertex.size();

			//dprint("\n\n\n");
			if (lp.region_vertex.size() > 1)
			{
				seed_vertex.push_back(ip.vl);
				dprint("seed vertex:", seed_index, ip.vl->v.idx(), ip.vl->id, ip.data);
			}
			else
			{
				//dprint("undrafted vertex:", ip.vl->v.idx(), ip.vl->id, ip.data);
			}
			if (lp.region_vertex.size() > 1)
			{
				for (auto rv : lp.region_vertex)
				{
					region_index[rv->id].push_back(seed_index);
					region_index[m4.another_layer(rv, 2)->id].push_back(seed_index);
				}
			}
			
			auto& region_face = lp.region_face;
			auto& apl = lp.all_pl;
			for (int i = 0; i < apl.size(); ++i)
			{
				if (apl[i].empty())
					continue;
				//all_plane_loop[i] = std::move(apl[i]);
			}

			for (auto rf : region_face)
			{
				constraint_flag[rf->f.idx()] = true;
				constraint_dir.col(rf->f.idx()) = lp.x_axis.col(rf->f.idx());
			}
#if 1
			if (lp.region_vertex.size() > 1)
			{
				cset.cylinders.push_back(cylinder());
				cylinder &cy = cset.cylinders.back();
				cy.id = cset.cylinders.size() - 1;
				cy.vertices = std::move(lp.region_vertex);
				cy.vertice_flag = std::move(lp.region_v_flag);
				cy.faces = std::move(lp.region_face);
				cy.face_flag = std::move(lp.region_f_flag);
			}
#endif
			if (ip.vl->v.idx() == 6440)
			{
				//break;
			}
		}
		//return;
		ProcessOverlap();
		WriteRegion(cset.cylinders, cf->crossfield, model_name);
		for (auto &cy : cset.cylinders)
		{
			for (auto rf : cy.faces)
			{
				constraint_flag[rf->f.idx()] = true;
				constraint_dir.col(rf->f.idx()) = cf->crossfield.col(rf->id);
			}
		}
		cf->setOuterConstraint(constraint_flag, constraint_dir);
		cf->setField();
		cf->initFieldInfo();
		cf->write_field();
		tr.out("repair field:");
		dprint("Construct Cylinder Done");
	}

	void LoopGen::OptimizeCylinder()
	{
		m4.set_base(mesh, cf);
		m4.init();
		m4.update();
		m4.set_weight();

		std::vector<std::vector<int>> vh_set; vh_set.reserve(cset.cylinders.size());
		std::vector<OpenMesh::Vec3d> dir; dir.reserve(cset.cylinders.size());
#if 1
		for (auto &cy : cset.cylinders)
		{
			int fid = cy.vertices.front()->hl->left;
			std::vector<int> vs; vs.reserve(cy.vertices.size());
			for (auto vl : cy.vertices)
				vs.push_back(vl->v.idx());
			vh_set.push_back(std::move(vs));
			dir.emplace_back(cf->crossfield(0, fid), cf->crossfield(1, fid), cf->crossfield(2, fid));
		}
#else
		ReadRegion(vh_set, dir, model_name);
#endif
		RecoverCylinder(vh_set, dir, false, false, false, false, false);

		auto &crossfield = cf->getCrossField();
		std::vector<std::vector<int>> region_index(m4.verticelayers.size());
		BoolVector constraint_flag(mesh->n_faces(), false);
		Eigen::Matrix3Xd constraint_dir(3, mesh->n_faces());
		seed_vertex.clear();
		for (auto &cy : cset.cylinders)
		{
			VertexLayer* vl = cy.vertices.front();
			LocalParametrization lp(m4, vl);
			//初始化初始区域，将读入的cylinder设置为new region
			lp.new_vertex = std::vector<VertexLayer*>(cy.vertices.begin() + 1, cy.vertices.end());
			lp.new_v_flag = cy.vertice_flag;
			lp.new_v_flag[vl->id] = false;
			lp.new_face = std::move(cy.faces);
			lp.new_f_flag = std::move(cy.face_flag);
			auto &x_axis = lp.x_axis;
			auto &y_axis = lp.y_axis;
			for (auto fl : lp.new_face)
			{
				int flid = fl->id;
				int fid = flid / 4;
				x_axis.col(fid) = crossfield.col(flid);
				y_axis.col(fid) = crossfield.col(fid * 4 + (flid + 1) % 4);
			}
			bool grow_flag[4] = { true, true, true, true };
			do
			{
				BoolVector visited_v = lp.region_v_flag;
				for (auto& ver : lp.new_vertex)
					visited_v[ver->id] = true;
				if (!ConstructRegionCut(vl, visited_v, lp.cut))
					break;
				lp.run(cf->getNormal());
			} while (SpreadSubRegion(lp, grow_flag));
			//dprint("vertices:", lp.region_vertex.size());
			if (lp.region_vertex.size() > 1)
				lp.modify_cut();

			int seed_index = seed_vertex.size();
			if (lp.region_vertex.size() > 1)
			{
				seed_vertex.push_back(vl);
				dprint("seed vertex:", seed_index, vl->v.idx(), vl->id);
			}
			else
			{
				//dprint("undrafted vertex:", ip.vl->v.idx(), ip.vl->id, ip.data);
			}
			if (lp.region_vertex.size() > 1)
			{
				for (auto rv : lp.region_vertex)
				{
					region_index[rv->id].push_back(seed_index);
					region_index[m4.another_layer(rv, 2)->id].push_back(seed_index);
				}
			}
			for (auto rf : lp.region_face)
			{
				constraint_flag[rf->f.idx()] = true;
				constraint_dir.col(rf->f.idx()) = lp.x_axis.col(rf->f.idx());
			}
#if 1
			if (lp.region_vertex.size() > 1)
			{
				cy.vertices = std::move(lp.region_vertex);
				cy.vertice_flag = std::move(lp.region_v_flag);
			}
#endif
		}
		ProcessOverlap();
		WriteRegion(cset.cylinders, cf->crossfield, model_name);
		cf->setOuterConstraint(constraint_flag, constraint_dir);
		cf->setField();
		cf->initFieldInfo();
		cf->write_field();
		tr.out("repair field:");
		dprint("Optimize Cylinder Done");
	}

	void LoopGen::IterateCylinder()
	{
		//参考现有cylinder进行顶点上loop的搜索
#if 1
		m4.set_base(mesh, cf);
		m4.init();
		m4.update();
		m4.set_weight();

		int nvl = m4.verticelayers.size();
		std::vector<std::vector<int>> vh_set; vh_set.reserve(cset.cylinders.size());
		std::vector<OpenMesh::Vec3d> dir; dir.reserve(cset.cylinders.size());
#if 0
		for (auto &cy : cset.cylinders)
		{
			int fid = cy.vertices.front()->hl->left;
			std::vector<int> vs; vs.reserve(cy.vertices.size());
			for (auto vl : cy.vertices)
				vs.push_back(vl->v.idx());
			vh_set.push_back(std::move(vs));
			dir.emplace_back(cf->crossfield(0, fid), cf->crossfield(1, fid), cf->crossfield(2, fid));
		}
#else
		InfoOnMesh.resize(m4.verticelayers.size());
		int count = 0;
		for (int i = 0; i < InfoOnMesh.size(); ++i)
		{
			InfoOnMesh[i].id = i;
			InfoOnMesh[i].energy = YYSS_INFINITE;
			switch (i % 4)
			{
			case 0:
				InfoOnMesh[i].dir = Forward;
				InfoOnMesh[i].plid = count;
				++count;
				break;
			case 1:
				InfoOnMesh[i].dir = Forward;
				InfoOnMesh[i].plid = count;
				break;
			case 2:
				InfoOnMesh[i].dir = Reverse;
				InfoOnMesh[i].plid = count - 1;
				break;
			case 3:
				InfoOnMesh[i].dir = Reverse;
				InfoOnMesh[i].plid = count;
				++count;
				break;
			}
		}
		pls.resize(nvl / 2);
		ReadRegion(vh_set, dir, model_name/* + "2"*/);
#endif
		RecoverCylinder(vh_set, dir, true, true, true, true, false);
#endif



		//return;
#if 1
		pls.resize(m4.verticelayers.size() / 2);
		ReadPlaneLoop(m4, pls, model_name + "2", mesh);
		similarity_energy.resize(m4.halfedgelayers.size(), YYSS_INFINITE);
		ReadEnergy(similarity_energy, model_name + "2");
#else
		std::vector<std::vector<VertexLayer*>> link_on_cylinder(nvl);
		double u = 0;
		double ru;
		auto adjust_u = [&]()
		{
			if (fabs(u - ru) > 0.5)
			{
				if (u > ru)
					u -= 1.0;
				else
					u += 1.0;
			}
		};
		for (auto &cy : cset.cylinders)
		{
			for (int i = 0; i < 2; ++i)
			{
				for (auto hl : cy.bounds[i])
				{
					VertexLayer* start_vl = m4.another_layer(&m4.verticelayers[hl->to], 1 + 2 * i);
					auto &link_p = link_on_cylinder[start_vl->id];
					link_p.push_back(start_vl);
					ru = cy.GetRegularU(hl->to);
					auto hl_begin = hl->oppo;
					auto hl_transfer = hl_begin;
					do
					{
						double u0 = cy.GetRegularU(hl_transfer->next->to);
						double u1 = cy.GetRegularU(hl_transfer->to);
						if (i)
							std::swap(u0, u1);
						if (u0 > u1 + 0.5)
						{
							if (ru > u0)
								u1 += 1.0;
							else if (ru < u1)
								u0 -= 1.0;
						}
						if (u0 < ru && u1 > ru)
						{
							hl_transfer = hl_transfer->next;
							break;
						}
						hl_transfer = hl_transfer->prev->oppo;
					} while (hl_transfer != hl_begin);

					link_p.push_back(m4.another_layer(&m4.verticelayers[hl_transfer->to], 1 + 2 * i));
					while (true)
					{
						hl_transfer = hl_transfer->oppo;
						u = cy.GetRegularU(hl_transfer->next->to);
						adjust_u();
						if ((u > ru) ^ i)
						{
							hl_transfer = hl_transfer->prev;
							VertexLayer* temp_vl = m4.another_layer(&m4.verticelayers[hl_transfer->to], 1 + 2 * i);
							if (temp_vl != link_p.back())
								link_p.push_back(temp_vl);
						}
						else
						{
							hl_transfer = hl_transfer->next;
							VertexLayer* temp_vl = m4.another_layer(&m4.verticelayers[hl_transfer->from], 1 + 2 * i);
							if (temp_vl != link_p.back())
								link_p.push_back(temp_vl);
						}
						if (cset.bound_edge_flag[hl_transfer->id / 8])
						{
							VertexLayer* temp_vl = m4.another_layer(&m4.verticelayers[hl_transfer->from], 1 + 2 * i);
							if (temp_vl != link_p.back())
								link_p.push_back(temp_vl);
							temp_vl = m4.another_layer(&m4.verticelayers[hl_transfer->to], 1 + 2 * i);
							if (temp_vl != link_p.back())
								link_p.push_back(temp_vl);
							break;
						}
					}
				}
			}
		}
		
		//2.构造cylinder外的点的loop
		for (auto &iov : InfoOnMesh)
			iov.energy = YYSS_INFINITE;
		int nv = mesh->n_vertices();
#pragma omp parallel for
		for (int i = 0; i < nv; ++i)
		{
			if (i != 1420)
			{
				//continue;
			}
			if (m4.sing_flag[i] || cset.has_vertex[i])
				continue;
			VertexLayer* vl = &m4.verticelayers[i * 4];
			std::vector<VertexLayer*> loop;
			if (CylinderBasedPLSearch(vl, loop, link_on_cylinder))
				RefineLoopByPlanarity(loop, pls[InfoOnMesh[vl->id].plid]);
			/*search_loop = std::move(loop);
			return;*/
			//dprint("size", loop.size());
			++vl;
			if (CylinderBasedPLSearch(vl, loop, link_on_cylinder))
				RefineLoopByPlanarity(loop, pls[InfoOnMesh[vl->id].plid]);
			//dprint("size", loop.size());
		}
		WritePlaneLoop(pls, model_name, mesh);

		InitializeSimilarityEnergy();
#endif
		for (auto &cy : cset.cylinders)
		{
			for (auto fl : cy.faces)
			{
				for (int j = 0; j < 8; ++j)
				{
					similarity_energy[fl->hl->id / 8 * 8 + j] = 0;
					similarity_energy[fl->hl->next->id / 8 * 8 + j] = 0;
					similarity_energy[fl->hl->prev->id / 8 * 8 + j] = 0;
				}
			}
			for (int i = 0; i < 2; ++i)
			{
				for (auto hl : cy.bounds[i])
				{
					HalfedgeLayer* hl_transfer = hl;
					do
					{
						for (int j = 0; j < 8; ++j)
						{
							similarity_energy[hl_transfer->id / 8 * 8 + j] = 0;
						}
						hl_transfer = hl_transfer->prev->oppo;
					} while (hl_transfer != hl);
				}
			}
		}
		for (auto &iov : InfoOnMesh)
			iov.energy = YYSS_INFINITE;
		for (auto &hl : m4.halfedgelayers)
		{
			InfoOnMesh[hl.from].energy = std::min(InfoOnMesh[hl.from].energy, similarity_energy[hl.id]);
			InfoOnMesh[hl.to].energy = std::min(InfoOnMesh[hl.to].energy, similarity_energy[hl.id]);
		}
		//return;
		WriteEnergy(similarity_energy, model_name);
		ConstructCylinder();
		WriteRegion(cset.cylinders, cf->crossfield, model_name);
		
	}

	bool LoopGen::CylinderBasedPLSearch(VertexLayer* vl, std::vector<VertexLayer*> &loop, std::vector<std::vector<VertexLayer*>> &link_on_cylinder)
	{
		layernode_pq pq;

		int nvl = m4.verticelayers.size();
		int vid;
		std::vector<double> distance(nvl, YYSS_INFINITE);
		std::vector<int> count(nvl, 0);
		std::vector<HalfedgeLayer*> prev(nvl, nullptr);
		BoolVector visited(nvl, false);

		HalfedgeLayer* hl_begin = vl->hl;
		HalfedgeLayer* hl_transfer = hl_begin;
		auto &crossfield = cf->getCrossField();
		do
		{
			double w = m4.weight(hl_transfer->id);
			if (w < YYSS_INFINITE)
			{
				vid = m4.verticelayers[hl_transfer->to].v.idx();
				if (!m4.sing_flag[vid])
				{
					int toid = hl_transfer->to;
					pq.emplace(toid, w, ++count[toid]);
					distance[toid] = w;
					prev[toid] = hl_transfer;
					//vid = m4.verticemap[vid];
					vid *= 4;
					for (int i = 0; i < 4; ++i)
					{
						if (vid + i == toid)
							continue;
						visited[vid + i] = true;
					}
				}
			}
			hl_transfer = hl_transfer->prev->oppo;
		} while (hl_transfer != hl_begin);

		BoolVector break_info(nvl, false);
		break_info[vl->id] = true;
		//VertexLayer* last_vl = nullptr;
		/*static int ess = 0;
		++ess;
		if (ess < 2)
		{
			int p = 0;
		}*/
		int cross_time = 0;
		while (true)
		{
			LayerNode ln;
			do
			{
				if (pq.empty())
					return false;
				ln = pq.top();
				pq.pop();
			} while (ln.count != count[ln.id]);
			int fromid = ln.id;
			//if (fromid == vl->id)
			if (break_info[fromid])
			{
				//last_vl = &m4.verticelayers[fromid];
				break;
			}

			if (link_on_cylinder[fromid].empty())
			{
				hl_begin = m4.verticelayers[fromid].hl;
				hl_transfer = hl_begin;
				do
				{
					double w = m4.weight(hl_transfer->id);
					if (w < YYSS_INFINITE)
					{
						vid = m4.verticelayers[hl_transfer->to].v.idx();
						if (!m4.sing_flag[vid] && !visited[hl_transfer->to])
						{
							int toid = hl_transfer->to;
							if (distance[fromid] + w < distance[toid])
							{
								distance[toid] = distance[fromid] + w;
								pq.emplace(toid, distance[toid], ++count[toid]);
								prev[toid] = hl_transfer;
								vid *= 4;
								for (int i = 0; i < 4; ++i)
								{
									if (vid + i == toid)
										continue;
									visited[vid + i] = true;
								}
							}
						}
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);
			}
			else
			{
				++cross_time;
				if (cross_time > 10 * cset.cylinders.size())
					return false;
				auto &link_path = link_on_cylinder[fromid];
				HalfedgeLayer* first_hl = m4.find_halfedge_layer(link_path[0], link_path[1]);
				for (int i = 0; i < 2; ++i)
				{
					VertexLayer* last_vl = link_path[link_path.size() - 1 - i];
					vid = last_vl->v.idx();
					int toid = last_vl->id;
					if (!visited[toid])
					{
						distance[toid] = distance[fromid];
						pq.emplace(toid, distance[toid], ++count[toid]);
						prev[toid] = first_hl;
						vid *= 4;
						for (int j = 0; j < 4; ++j)
						{
							if (vid + j == toid)
								continue;
							visited[vid + i] = true;
						}
					}
				}
			}
		}
		/*std::ofstream file_writer;
		file_writer.open("C://Users//123//Desktop//vscode compare//old.txt");
		if (file_writer.fail()) {
			return false;
		}*/
		loop.clear();
		loop.push_back(vl);
		VertexLayer* prev_layer = &m4.verticelayers[prev[vl->id]->from];
		while (prev_layer != vl)
		{
			loop.push_back(prev_layer);
			//file_writer << prev_layer->id << std::endl;
			if (prev[prev_layer->id]->to != prev_layer->id)
			{
				auto &link_path = link_on_cylinder[prev[prev_layer->id]->from];
				for (int i = link_path.size() - 2; i > 0; --i)
					loop.push_back(link_path[i]);
			}
			prev_layer = &m4.verticelayers[prev[prev_layer->id]->from];
			if (loop.size() > 1.0e7)
				break;
		}
		loop.push_back(vl);
		//file_writer.close();
		std::reverse(loop.begin(), loop.end());
		return true;
	}

	void LoopGen::RecoverCylinder(std::vector<std::vector<int>> &vh_set, std::vector<OpenMesh::Vec3d> &dir,
		bool set_cut, bool set_bound, bool set_parameter, bool set_flag, bool set_intersection)
	{
		int size = vh_set.size();
		cset.cylinders.swap(std::vector<cylinder>());
		cset.cylinders.resize(size);
		double halfSquare2 = sqrt(2)*0.5;
		//恢复柱体数据
		for (int i = 0; i < size; ++i)
		{
			BoolVector vh_flag(mesh->n_vertices(), false);
			for (auto vh : vh_set[i])
				vh_flag[vh] = true;
			//VertexLayer* seed_vl = &m4.verticelayers[m4.verticemap[vh_set[i].front()]];
			VertexLayer* seed_vl = &m4.verticelayers[vh_set[i].front() * 4];
			auto &vec = m4.cf->crossfield.col(seed_vl->hl->left);
			double dot_ = dir[i][0] * vec(0) + dir[i][1] * vec(1) + dir[i][2] * vec(2);
			if (dot_ > -halfSquare2 && dot_ < halfSquare2)
				seed_vl = m4.another_layer(seed_vl, 1);

			//恢复顶点集合
			int begin_ = 0;
			int end_ = 1;
			auto &vls = cset.cylinders[i].vertices;
			auto &vl_flag = cset.cylinders[i].vertice_flag;
			vls.push_back(seed_vl);
			vl_flag.resize(m4.verticelayers.size(), false);
			vl_flag[seed_vl->id] = true;
			do
			{
				for (int j = begin_; j < end_; ++j)
				{
					HalfedgeLayer* hl_begin = vls[j]->hl;
					HalfedgeLayer* hl_transfer = hl_begin;
					do
					{
						VertexLayer* vl = &m4.verticelayers[hl_transfer->to];
						if (vh_flag[vl->v.idx()] && !vl_flag[vl->id] && !m4.sing_flag[vl->v.idx()])
						{
							vls.push_back(vl);
							vl_flag[vl->id] = true;
						}
						hl_transfer = hl_transfer->prev->oppo;
					} while (hl_transfer != hl_begin);
				}
				begin_ = end_;
				end_ = vls.size();
			} while (begin_ != end_);
			cset.cylinders[i].id = i;

			//恢复面集合
			auto &fls = cset.cylinders[i].faces;
			auto &fl_flag = cset.cylinders[i].face_flag;
			fl_flag.resize(m4.facelayers.size(), false);
			for (auto vl : vls)
			{
				HalfedgeLayer* hl_begin = vl->hl;
				HalfedgeLayer* hl_transfer = hl_begin;
				do
				{
					if (!fl_flag[hl_transfer->left])
					{
						if (vl_flag[hl_transfer->to] && vl_flag[hl_transfer->next->to])
						{
							fls.push_back(&m4.facelayers[hl_transfer->left]);
							fl_flag[hl_transfer->left] = true;
						}
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);
			}

			//构造cut
			if(set_cut)
				ConstructRegionCut(seed_vl, vl_flag, cset.cylinders[i].cut);

			//初始化bound
			if(set_bound)
				cset.cylinders[i].set_bound();

			//计算参数化
			//这个参数化代码好像是有bug
			if(set_parameter)
				cset.cylinders[i].parametrize(m4, cf->getNormal());
		}
		
		if (set_flag)
		{
			cset.bound_edge_flag.resize(mesh->n_edges(), false);
			cset.vertex_bound_index.resize(mesh->n_vertices(), std::pair<int, int>(-1, -1));
			cset.has_vertex.resize(mesh->n_vertices(), false);
			cset.has_face.resize(mesh->n_faces(), false);
			//恢复柱体数据
			for (auto &cy : cset.cylinders)
			{
				for (auto vl : cy.vertices)
					cset.has_vertex[vl->v.idx()] = true;
				for (auto fl : cy.faces)
					cset.has_face[fl->f.idx()] = true;
				for (int j = 0; j < 2; ++j)
				{
					for (auto hl : cy.bounds[j])
					{
						cset.bound_edge_flag[hl->id / 8] = true;
						cset.vertex_bound_index[hl->to / 4].first = cy.id;
						cset.vertex_bound_index[hl->to / 4].second = j;
					}
				}
			}
		}
		if (set_intersection)
		{
			std::vector<std::vector<int>> intersection(mesh->n_vertices());
			for (auto &cy : cset.cylinders)
			{
				for (auto vl : cy.vertices)
				{
					intersection[vl->v.idx()].push_back(cy.id);
				}
			}
			auto &tangential_intersection = cset.tangential_intersection;
			tangential_intersection.resize(cset.cylinders.size());
			for (auto &ti : tangential_intersection)
				ti.resize(cset.cylinders.size(), false);
			for (int k = 0; k < intersection.size(); ++k)
			{
				int size_ = intersection[k].size();
				if (size_ < 2)
					continue;
				if (size_ > 2)
					dprint("\nthere is a vertex in three different region\n");
				for (int i = 0; i < size_; ++i)
				{
					for (int j = i + 1; j < size_; ++j)
					{
						tangential_intersection[intersection[k][i]][intersection[k][j]] = true;
						tangential_intersection[intersection[k][j]][intersection[k][i]] = true;
					}
				}
			}
		}
	}

	bool LoopGen::RefineLoopByParametrization(VertexLayer* vl, LocalParametrization& lp, BoolVector& visited_v, BoolVector& visited_f)
	{
		double v_para = lp.GetV(vl->id);
		auto hl_begin = vl->hl;
		auto hl_transfer = hl_begin;
		if (!visited_f[hl_transfer->left])
			return false;
		double s0 = lp.GetV(hl_transfer->to) - v_para;
		double distance[2];
		PointOnHalfedgeLayer pohl[2];
		int id[2] = { 0,0 };

		double s1 = lp.GetV(hl_transfer->next->to) - v_para;
		do
		{
			if (visited_v[hl_transfer->next->to])
				s1 = lp.GetV(hl_transfer->next->to) - v_para;
			if (visited_f[hl_transfer->left])
			{
				if (s0*s1 < 0)
				{
					if (s0 < 0)
					{
						pohl[0].hl = hl_transfer->next;
						pohl[0].c = s1 / (s1 - s0);
						distance[0] = s0; distance[1] = s1;
						++id[0];
					}
					else
					{
						pohl[1].hl = hl_transfer->next->oppo;
						pohl[1].c = s0 / (s0 - s1);
						++id[1];
					}
				}
			}
			s0 = s1;
			hl_transfer = hl_transfer->prev->oppo;
		} while (hl_transfer != hl_begin);

		if (id[0] != 1 || id[1] != 1)
			return false;

		PlaneLoop planar_loop;
		planar_loop.push_back(pohl[0]);
		auto hl = pohl[0].hl;

		while (hl != pohl[1].hl)
		{
			int vlid = hl->oppo->next->to;
			if (!visited_v[vlid])
				return false;
			double s = lp.GetV(vlid) - v_para;
			if (s > 0)
			{
				hl = hl->oppo->next;
				distance[1] = s;
			}
			else
			{
				hl = hl->oppo->prev;
				distance[0] = s;
			}
			planar_loop.emplace_back(hl, distance[1] / (distance[1] - distance[0]));
		}

		lp.all_pl[vl->id] = std::move(planar_loop);
		return true;
	}

	double LoopGen::LoopLenGrad(std::vector<VertexLayer*>& vertex_set, LocalParametrization& lp, BoolVector& vertex_flag, int growDir)
	{
		int vsn = vertex_set.size();
		std::vector<double> len(vsn, -1.0);
		auto& grow_dir = lp.grow_dir;
#pragma omp parallel for
		for (int i = 0; i < vsn; ++i)
		{
			if (grow_dir[vertex_set[i]->id] != growDir)
				continue;
			const auto& pl = lp.all_pl[vertex_set[i]->id];
			int pln = pl.size();
			Vec3d pos[2] = { mesh->point(vertex_set[i]->v),
				pl.front().c * mesh->point(m4.verticelayers[pl.front().hl->from].v) +
				(1 - pl.front().c) * mesh->point(m4.verticelayers[pl.front().hl->to].v) };
			len[i] = 0;
			len[i] += (pos[0] - pos[1]).norm();
			for (int j = 1; j < pln; ++j)
			{
				pos[(j + 1) & 1] = pl[j].c * mesh->point(m4.verticelayers[pl[j].hl->from].v)
					+ (1 - pl[j].c) * mesh->point(m4.verticelayers[pl[j].hl->to].v);
				len[i] += (pos[0] - pos[1]).norm();
			}
			len[i] += (mesh->point(vertex_set[i]->v) - pl.back().c * mesh->point(m4.verticelayers[pl.back().hl->from].v)
				- (1 - pl.back().c) * mesh->point(m4.verticelayers[pl.back().hl->to].v)).norm();
		}
		double minmax_len[2] = { YYSS_INFINITE, 0.0 };
		int minmax_id[2] = { -1, -1 };
		for (int i = 0; i < vsn; ++i)
		{
			if (len[i] < 0.0)
				continue;
			if (len[i] < minmax_len[0])
			{
				minmax_len[0] = len[i];
				minmax_id[0] = vertex_set[i]->id;
			}
			if (len[i] > minmax_len[1])
			{
				minmax_len[1] = len[i];
				minmax_id[1] = vertex_set[i]->id;
			}
		}
		if (minmax_id[0] < 0 || minmax_id[1] < 0)
			return YYSS_INFINITE;
		if (minmax_id[0] == minmax_id[1])
			return YYSS_INFINITE;
		//计算 minmax_id 中两条 loop 在 u = 0.5 处的点的距离
		Vec3d pos[2];
		auto assembleU = [&](const PointOnHalfedgeLayer& poh)
		{
			double t0 = lp.GetRegularU(poh.hl->from);
			double t1 = lp.GetRegularU(poh.hl->to);
			if (fabs(t0 - t1) > 0.5) { if (t0 < t1) { t0 += 1.0; } else { t1 += 1.0; } }
			double u = poh.c * t0 + (1 - poh.c) * t1;
			return u - std::floor(u);
		};
		auto setdata = [&](int i)
		{
			double up[2];
			const auto& pl = lp.all_pl[minmax_id[i]];
			int pln = pl.size();
			up[0] = assembleU(pl.front());
			for (int j = 1; j < pln; ++j)
			{
				up[j & 1] = assembleU(pl[j]);
				if ((up[0] - 0.5) * (up[1] - 0.5) <= 0.0 && fabs(up[0] - up[1]) < 0.5)
				{
					Vec3d p0 = pl[j - 1].c * mesh->point(m4.verticelayers[pl[j - 1].hl->from].v)
						+ (1 - pl[j - 1].c) * mesh->point(m4.verticelayers[pl[j - 1].hl->to].v);
					Vec3d p1 = pl[j].c * mesh->point(m4.verticelayers[pl[j].hl->from].v)
						+ (1 - pl[j].c) * mesh->point(m4.verticelayers[pl[j].hl->to].v);
					double lambda = (0.5 - up[(j + 1) & 1]) / fabs(up[0] - up[1]);
					pos[i] = (1 - lambda) * p0 + lambda * p1;
					return true;
				}
			}
			return false;
		};
		setdata(0);
		setdata(1);
		if ((pos[0] - pos[1]).norm() <= YYSS_FAIRLY_SMALL)
			return YYSS_INFINITE;
		return (minmax_len[1] - minmax_len[0]) / (pos[0] - pos[1]).norm();
	}

	void LoopGen::AssembleIOVLoopEnergy(vl_pair_pq& pq)
	{
		for (int i = 0; i < mesh->n_vertices(); ++i)
		{
			if (m4.sing_flag[i])
				continue;
			if (cset.has_vertex.size() != 0 && cset.has_vertex[i])
				continue;
			for (int j = 0; j < 2; ++j)
			{
				VertexLayer* vl = &m4.verticelayers[i * 4 + j];
				auto& iov = InfoOnMesh[vl->id];
				auto& pl = pls[iov.plid];
				if (pl.empty())
					continue;
				double sum = iov.energy;
				for (auto pohl : pl)
				{
					sum += pohl.c * InfoOnMesh[pohl.hl->from].energy + (1 - pohl.c) * InfoOnMesh[pohl.hl->to].energy;
				}
				pq.emplace(vl, sum / (1 + pl.size()));
			}
		}
#if PRINT_DEBUG_INFO
		dprint("计算圈上的平均能量");
#endif
	}

	bool LoopGen::CheckCylinderTopology(std::vector<VertexLayer*>& vertex_set, BoolVector& vs_flag, std::vector<int>& grow_dir)
	{
		//check if empty
		if (vertex_set.empty())
			return false;
		//check connectivity
		BoolVector visited(m4.verticelayers.size(), false);
		std::queue<VertexLayer*> tree;
		tree.push(vertex_set.front());
		visited[vertex_set.front()->id] = true;
		int count = 0;
		while (!tree.empty())
		{
			auto vl = tree.front();
			tree.pop();
			++count;
			auto hl_begin = vl->hl;
			auto hl_transfer = hl_begin;
			do
			{
				if (vs_flag[hl_transfer->to] && !visited[hl_transfer->to])
				{
					tree.push(&m4.verticelayers[hl_transfer->to]);
					visited[hl_transfer->to] = true;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}
		if (count < vertex_set.size())
			return false;

		//check manifold
		BoolVector fs_flag(m4.facelayers.size(), false);
		BoolVector bv_flag(m4.verticelayers.size(), false);
		for (const auto& vl : vertex_set)
		{
			auto hl_begin = vl->hl;
			auto hl_transfer = hl_begin;
			do
			{
				if (!fs_flag[hl_transfer->left])
				{
					bool flag = true;
					HalfedgeLayer* hb = m4.facelayers[hl_transfer->left].hl;
					auto ht = hb;
					do
					{
						if (!vs_flag[ht->from])
						{
							flag = false;
							break;
						}
						ht = ht->next;
					} while (ht != hb);
					if (flag)
						fs_flag[hl_transfer->left] = true;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);

			if (bv_flag[vl->id])
				continue;
			hl_begin = vl->hl;
			hl_transfer = hl_begin;
			do
			{
				if (!vs_flag[hl_transfer->to])
				{
					bv_flag[vl->id] = true;
					break;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}
		for (const auto& vl : vertex_set)
		{
			auto hl_begin = vl->hl;
			auto hl_transfer = hl_begin;
			do
			{
				if (vs_flag[hl_transfer->to])
				{
					HalfedgeLayer* hl_ = m4.find_halfedge_layer(vl, &m4.verticelayers[hl_transfer->to]);
					if (!fs_flag[hl_->left] && !fs_flag[hl_->oppo->left])
						return false;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}

		//check boundary number and set grow_dir in passing
		BoolVector be_flag(m4.halfedgelayers.size() / 2, false);
		count = 0;
		for (const auto& vl : vertex_set)
		{
			if (!bv_flag[vl->id])
				continue;
			auto he = vl->hl;
			while (!fs_flag[he->left] || fs_flag[he->oppo->left])
				he = he->prev->oppo;
			if (be_flag[he->id / 2])
				continue;
			HalfedgeLayer* hl_transfer = he;
			do
			{
				grow_dir[hl_transfer->to] = count;
				be_flag[hl_transfer->id / 2] = true;
				hl_transfer = hl_transfer->oppo;
				while (!fs_flag[hl_transfer->left] || fs_flag[hl_transfer->oppo->left])
					hl_transfer = hl_transfer->prev->oppo;
			} while (he != hl_transfer);
			++count;
			if (count > 2)
				return false;
		}
		if (count < 2)
			return false;
		for (const auto& vl : vertex_set)
		{
			if (grow_dir[vl->id] == -1)
				grow_dir[vl->id] = 0;
		}
		return true;
	}

	void LoopGen::ProcessOverlap()
	{
		std::vector<std::vector<int>> region_index(m4.verticelayers.size());
		for (auto &cy : cset.cylinders)
		{
			for (auto vl : cy.vertices)
			{
				region_index[vl->id].push_back(cy.id);
				region_index[m4.another_layer(vl, 2)->id].push_back(cy.id);
			}
		}
		if (cset.cylinders.size() < 2)
			return;
		int nvl = region_index.size();
		for (int i = 0; i < nvl; ++i)
		{
			int ris = region_index[i].size();
			while (ris > 1)
			{
				int parent = region_index[i][ris - 2];
				int child = region_index[i][ris - 1];
				auto &pc = cset.cylinders[parent];
				auto &cc = cset.cylinders[child];
				if (pc.vertice_flag[i] && cc.vertice_flag[i])
				{
					for (auto vl : cc.vertices)
					{
						if (pc.vertice_flag[vl->id])
							continue;
						pc.vertices.push_back(vl);
						pc.vertice_flag[vl->id] = true;
					}
				}
				else
				{
					VertexLayer* oppo_vl = nullptr;
					for (auto vl : cc.vertices)
					{
						oppo_vl = m4.another_layer(vl, 2);
						if (pc.vertice_flag[oppo_vl->id])
							continue;
						pc.vertices.push_back(oppo_vl);
						pc.vertice_flag[oppo_vl->id] = true;
					}
				}
				cc.id = -1;
				for (int j = 0; j < nvl; ++j)
				{
					if (region_index[j].empty())
						continue;
					if (region_index[j].back() == child)
						region_index[j].pop_back();
				}
				ris = region_index[i].size();
			}
		}

		std::vector<cylinder> new_cylinders;
		for (auto &cy : cset.cylinders)
		{
			if (cy.id == -1)
				continue;
			new_cylinders.push_back(std::move(cy));
			//new_cylinders.push_back(cy);
			new_cylinders.back().id = new_cylinders.size() - 1;
		}
		cset.cylinders = std::move(new_cylinders);
	}

	void LoopGen::ReLoop()
	{
		m4.set_base(mesh, cf);
		m4.init();
		m4.update();
		m4.set_weight();

		std::vector<std::vector<int>> vh_set;
		std::vector<OpenMesh::Vec3d> dir;
		ReadRegion(vh_set, dir, model_name);
		RecoverCylinder(vh_set, dir, true, true, false, true, true);

		//从柱体边界出发搜索路径
		//1.用平面性来决定路径搜索方向
		//2.用场来决定路径搜索方向
#if 0
		int nv = mesh->n_vertices() / 10;
		const auto &normal = cf->getNormal();
		const auto &crossfield = cf->getCrossField();
		cset.all_path.resize(m4.verticelayers.size());
		Eigen::Vector4d plane; plane.setZero();
		double fromdis = 0;
		double todis = 0;
		auto sign_dist = [&](OpenMesh::Vec3d &p)
		{
			return plane(0)*p[0] + plane(1)*p[1] + plane(2)*p[2] + plane(3);
		};

		//vl_pair_pq path_pq;
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				auto &bound = cset.cylinders[i].bounds[j];
				int shift = j == 0 ? 3 : 1;
				for (auto hl : bound)
				{
					VertexLayer* vl = m4.another_layer(&m4.verticelayers[hl->to], shift);
					HalfedgeLayer* hl_begin = vl->hl;
					HalfedgeLayer* hl_transfer = hl_begin;
					Eigen::Vector3d field_dir; field_dir.setZero();
					Eigen::Vector3d normal_dir; normal_dir.setZero();
					do
					{
						if (cset.cylinders[i].face_flag[m4.another_layer(hl_transfer, 4 - shift)->left])
						{
							field_dir += crossfield.col(hl_transfer->left);
							normal_dir += normal.col(hl_transfer->left / 4);
						}
						hl_transfer = hl_transfer->prev->oppo;
					} while (hl_transfer != hl_begin);
					plane.head(3) = field_dir.cross(normal_dir).normalized();
					auto start_pos = mesh->point(vl->v);
					plane(3) = -(plane(0)*start_pos[0] + plane(1)*start_pos[1] + plane(2)*start_pos[2]);
					/*do
					{
						plane.head(3) = crossfield.col(hl_transfer->left).cross(normal.col(hl_transfer->left / 4));
						auto &pos = mesh->point(vl->v);
						plane(3) = -(pos[0] * plane(0) + pos[1] * plane(1) + pos[2] * plane(2));
						auto &frompos = mesh->point(m4.verticelayers[hl_transfer->to].v);
						auto &topos = mesh->point(m4.verticelayers[hl_transfer->next->to].v);
						fromdis = sign_dist(frompos);
						todis = sign_dist(topos);
						if (fromdis > 0 && todis < 0)
						{
							hl_transfer = hl_transfer->next;
							break;
						}
						hl_transfer = hl_transfer->prev->oppo;
					} while (hl_transfer != hl_begin);*/
					auto &path = cset.all_path[vl->id];//注意此顶点已经不是region上的顶点了
					do
					{
						if (sign_dist(mesh->point(m4.verticelayers[hl_transfer->to].v)) >= 0
							&& sign_dist(mesh->point(m4.verticelayers[hl_transfer->next->to].v)) <= 0)
						{
							hl_transfer = hl_transfer->next;
							break;
						}
						hl_transfer = hl_transfer->prev->oppo;
					} while (hl_transfer != hl_begin);
					path.emplace_back(hl_transfer, todis / (todis - fromdis));
					int itertimes = 0;
					while (true)
					{
						hl_transfer = hl_transfer->oppo;
						auto pos = path.back().point(m4);
						double pos_dis = sign_dist(mesh->point(m4.verticelayers[hl_transfer->next->to].v));
						if (pos_dis >= 0)
						{
							fromdis = pos_dis;
							hl_transfer = hl_transfer->prev;
							path.emplace_back(hl_transfer, todis / (todis - fromdis));
						}
						else
						{
							todis = pos_dis;
							hl_transfer = hl_transfer->next;
							path.emplace_back(hl_transfer, todis / (todis - fromdis));
						}
						if (cset.bound_edge_flag[hl_transfer->id / 8])
							break;
						if (++itertimes > nv)
						{
							path.clear(); break;
						}
					}
					if (path.empty())
						continue;
					auto &fromb = cset.vertex_bound_index[vl->v.idx()];
					auto &tob = cset.vertex_bound_index[path.back().hl->to / 4];
					if (fromb == tob || cset.tangential_intersection[fromb.first][tob.first])
						path.clear();
				}
			}
		}
#else
		const auto &normal = cf->getNormal();
		const auto &crossfield = cf->getCrossField();
		cset.all_path.resize(m4.verticelayers.size());
		Eigen::Vector4d plane;
		double fromdis = 0;
		double todis = 0;
		auto sign_dist = [&](OpenMesh::Vec3d &p)
		{
			return plane(0)*p[0] + plane(1)*p[1] + plane(2)*p[2] + plane(3);
		};

		//vl_pair_pq path_pq;
		int size = vh_set.size();
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				auto &bound = cset.cylinders[i].bounds[j];
				int shift = j == 0 ? 3 : 1;
				for (auto hl : bound)
				{
					VertexLayer* vl = m4.another_layer(&m4.verticelayers[hl->to], shift);
					HalfedgeLayer* hl_begin = vl->hl;
					HalfedgeLayer* hl_transfer = hl_begin;
					do
					{
						plane.head(3) = crossfield.col(hl_transfer->left).cross(normal.col(hl_transfer->left / 4));
						auto &pos = mesh->point(vl->v);
						plane(3) = -(pos[0] * plane(0) + pos[1] * plane(1) + pos[2] * plane(2));
						auto &frompos = mesh->point(m4.verticelayers[hl_transfer->to].v);
						auto &topos = mesh->point(m4.verticelayers[hl_transfer->next->to].v);
						fromdis = sign_dist(frompos);
						todis = sign_dist(topos);
						if (fromdis > 0 && todis < 0)
						{
							hl_transfer = hl_transfer->next;
							break;
						}
						hl_transfer = hl_transfer->prev->oppo;
					} while (hl_transfer != hl_begin);
					auto &path = cset.all_path[vl->id];//注意此顶点已经不是region上的顶点了
					if (cset.has_vertex[hl_transfer->from / 4] && cset.has_vertex[hl_transfer->to / 4])
						continue;
					path.emplace_back(hl_transfer, todis / (todis - fromdis));

					int itertimes = 0;
					auto leftpos = mesh->point(vl->v);
					double l = 0;
					while (!cset.bound_edge_flag[hl_transfer->id / 8])
					{
						hl_transfer = hl_transfer->oppo;
						plane.head(3) = crossfield.col(hl_transfer->left).cross(normal.col(hl_transfer->left / 4));
						auto pos = path.back().point(m4);
						plane(3) = -(pos[0] * plane(0) + pos[1] * plane(1) + pos[2] * plane(2));
						double dis = sign_dist(mesh->point(m4.verticelayers[hl_transfer->next->to].v));
						if (dis > 0)
						{
							fromdis = dis;
							hl_transfer = hl_transfer->prev;
							todis = sign_dist(mesh->point(m4.verticelayers[hl_transfer->to].v));
							path.emplace_back(hl_transfer, todis / (todis - fromdis));
						}
						else
						{
							todis = dis;
							hl_transfer = hl_transfer->next;
							fromdis = sign_dist(mesh->point(m4.verticelayers[hl_transfer->from].v));
							path.emplace_back(hl_transfer, todis / (todis - fromdis));
						}
						l += (leftpos - pos).norm();
						leftpos = pos;
						if (++itertimes > 1000)
						{
							path.clear();
							break;
						}
					}

					if (path.empty())
						continue;
					auto &fromb = cset.vertex_bound_index[vl->v.idx()];
					auto &tob = cset.vertex_bound_index[path.back().hl->to / 4];
					if (fromb == tob || cset.tangential_intersection[fromb.first][tob.first])
						path.clear();
				}
			}
		}
#endif
		for (auto &cy : cset.cylinders)
		{
			cy.handle_to_layer.resize(mesh->n_vertices(), std::pair<int, int>(-1, -1));
			for (int i = 0; i < 2; ++i)
			{
				int shift = i == 0 ? 3 : 1;
				for (auto hl : cy.bounds[i])
				{
					cy.handle_to_layer[hl->to / 4].first = hl->to;
					cy.handle_to_layer[hl->to / 4].second = m4.another_layer(&m4.verticelayers[hl->to], shift)->id;
				}
			}
		}

		//OptimizeDisk(path_pq);
		OptimizeDisk();
	}

	void LoopGen::OptimizeDisk()
	{
		auto assemblePath = [&](PlaneLoop &pl, OpenMesh::Vec3d &start, Eigen::Matrix3Xd &path)
		{
			path.resize(3, 1 + pl.size());
			path.col(0) << start[0], start[1], start[2];
			int c = 0;
			for (auto itr = pl.begin(); itr != pl.end(); ++itr)
			{
				auto pos = itr->point(m4);
				path.col(++c) << pos[0], pos[1], pos[2];
			}
		};
		int cy_size = cset.cylinders.size();
		vl_pair_pq pq;
		//std::vector<double> path_similarity_energy(mesh->n_edges(), YYSS_INFINITE);
#if 1
		const auto &crossfield = cf->getCrossField();
		for (int i = 0; i < cy_size; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				int shift = j == 0 ? 3 : 1;
				auto &bound = cset.cylinders[i].bounds[j];
				for (int k = 0; k < bound.size(); ++k)
				{
					VertexLayer* vl = m4.another_layer(&m4.verticelayers[bound[k]->from], shift);
					auto &path = cset.all_path[vl->id];
					if (path.empty())
						continue;
					Eigen::Vector3d from_dir; from_dir.setZero();
					Eigen::Vector3d to_dir; to_dir.setZero();
					for (auto hl = bound[k]->oppo; hl != (k == 0 ? bound.back() : bound[k - 1]); hl = hl->next->oppo)
					{
						from_dir += crossfield.col(hl->left);
					}
					from_dir.normalize();
					to_dir = crossfield.col(m4.another_layer(&m4.facelayers[path.back().hl->oppo->left], shift)->id);

					double e = 0;
					double length = 0;
					Eigen::Matrix3Xd pos;
					assemblePath(path, mesh->point(mesh->vertex_handle(bound[k]->from / 4)), pos);
					for (int m = 0; m < pos.cols() - 1; ++m)
					{
						auto dif = pos.col(m + 1) - pos.col(m);
						e += fabs(dif.dot(from_dir)) + fabs(dif.dot(to_dir));
						length += dif.norm();
					}
					pq.emplace(m4.another_layer(&m4.verticelayers[bound[k]->from], shift), e * 0.5 / length);
				}
			}
		}
//#else
		edge_path_similarity_energy.resize(mesh->n_edges(), YYSS_INFINITE);
		vert_path_similarity_energy.resize(mesh->n_vertices(), YYSS_INFINITE);
		all_vertice_path.resize(mesh->n_vertices());
		for (int i = 0; i < cy_size; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				/*if (i == 2 && j == 1)
				{
					int p = 11745;
				}*/
				int shift = j == 0 ? 3 : 1;
				auto &bound = cset.cylinders[i].bounds[j];
				Eigen::Matrix3Xd path0, path1;
				Eigen::Matrix3Xd* pointer0 = &path0;
				Eigen::Matrix3Xd* pointer1 = &path1;
				VertexLayer* vl = m4.another_layer(&m4.verticelayers[bound[0]->from], shift);
				assemblePath(cset.all_path[vl->id], mesh->point(mesh->vertex_handle(bound[0]->from / 4)), path0);
				//std::vector<double> path_similarity_energy(bound.size(), YYSS_INFINITE);
				for (int k = 0; k < bound.size(); ++k)
				{
					/*if (k == 16)
					{
						int p = 0;
					}*/
					vl = m4.another_layer(&m4.verticelayers[bound[k]->to], shift);
					assemblePath(cset.all_path[vl->id], mesh->point(mesh->vertex_handle(bound[k]->to / 4)), *pointer1);
					Eigen::Matrix3Xd* ptr = pointer0; pointer0 = pointer1; pointer1 = ptr;
					VertexLayer* vl_another = m4.another_layer(&m4.verticelayers[bound[k]->from], shift);
					if (cset.all_path[vl->id].empty() || cset.all_path[vl_another->id].empty()
						|| cset.vertex_bound_index[cset.all_path[vl->id].back().hl->to / 4] !=
						cset.vertex_bound_index[cset.all_path[vl_another->id].back().hl->to / 4])
						edge_path_similarity_energy[bound[k]->id / 8] = YYSS_INFINITE;
					else
						edge_path_similarity_energy[bound[k]->id / 8] = EvaluatePathSimilarity(*pointer0, *pointer1);
				}
				for (int k = 0; k < bound.size(); ++k)
				{
					/*if (bound[k]->to/4 == 11745)
					{
						int p = 11745;
					}*/
					double e = 0.5*(edge_path_similarity_energy[bound[k]->id / 8] + 
						edge_path_similarity_energy[bound[(k + 1) % bound.size()]->id / 8]);
					vert_path_similarity_energy[bound[k]->to / 4] = e;
					//if (e < YYSS_INFINITE)
					//	pq.emplace(m4.conj_vl(&m4.verticelayers[bound[k]->to], shift), e);
				}
			}
		}
		similarity_energy = edge_path_similarity_energy;
#endif


		//return;


		BoolVector set_flag(mesh->n_vertices(), false);
		old_face_flag.resize(mesh->n_faces(), false);
		while (true)
		{
			vl_pair vp;
			vp = pq.top(); pq.pop();
			while (set_flag[vp.vl->id / 4])
			{
				if (pq.empty())
				{
					vp.data = YYSS_INFINITE;
					break;
				}
				vp = pq.top(); pq.pop();
			}
			//dprint(vp.data);
			if (vp.data > 0.1)
				break;


			bool grow_flag[2] = { true,true };
			dset.disks.push_back(disk());
			auto &dk = dset.disks.back();
			dk.from_bound = cset.vertex_bound_index[vp.vl->v.idx()];
			dk.to_bound = cset.vertex_bound_index[cset.all_path[vp.vl->id].back().hl->to / 4];
			dprint(vp.vl->id, vp.vl->v.idx(), dk.from_bound.first, dk.from_bound.second, dk.to_bound.first, dk.to_bound.second);
			/*if (vp.vl->id == 142502)
				break;*/
			//seed_vertex.push_back(vp.vl);
			if (vp.vl->id == 23937)
			{
				int p = 0;
			}

			ConstructInitialRegion(vp.vl, dk);

			if (vp.vl->id == 89113)
			{
				int p = 0;
				//OptimizeDiskField(dk);
			}
			//else
			{
				do
				{
					OptimizeDiskField(dk);
				} while (SpreadSubRegion(dk, grow_flag));
			}
			if (!dk.region_face.empty())
			{
				dprint(vp.vl->id, vp.vl->v.idx(), vp.data);
				seed_vertex.push_back(vp.vl);
			}

			for (auto rv : dk.region_vertex)
			{
				set_flag[rv.idx()] = true;
			}
			new_vert_flag = dk.new_v_flag;
			old_vert_flag = dk.region_v_flag;
			new_face_flag = dk.new_f_flag;
#if 0
			old_face_flag = dk.region_f_flag;

#else
			for (auto rf : dk.region_face)
				old_face_flag[rf.idx()] = true;
#endif

			xaxis = dk.x_axis;


			if (!dk.region_face.empty())
			{
				for (auto v : dk.region_vertex)
				{
					for (auto voh : mesh->voh_range(v))
					{
						if (dk.region_f_flag[voh.face().idx()] ^ dk.region_f_flag[voh.opp().face().idx()])
							dk.bounds.push_back(voh);
					}
				}
			}
			dk.initial_path = all_vertice_path[dk.region_vertex.front().idx()];
			if (dk.region_face.empty())
				dset.disks.pop_back();
			if (vp.vl->id == 89113)
			{
				//break;
			}
		}
		BoolVector constraint_flag(mesh->n_faces(), false);
		Eigen::Matrix3Xd constraint_dir(3, mesh->n_faces());
		for (auto &dk : dset.disks)
		{
			for (auto fh : dk.region_face)
			{
				constraint_flag[fh.idx()] = true;
				constraint_dir.col(fh.idx()) = dk.x_axis.col(fh.idx());
			}
		}
		for (auto &cy : cset.cylinders)
		{
			for (auto fh : cy.faces)
			{
				constraint_flag[fh->f.idx()] = true;
				constraint_dir.col(fh->f.idx()) = cf->crossfield.col(fh->id);
			}
		}
		//cf->setOuterConstraint(constraint_flag, constraint_dir);
		//cf->setField();
		//cf->initFieldInfo();
	}

	void LoopGen::ConstructInitialRegion(VertexLayer* vl, disk &dk)
	{
		auto &pl = cset.all_path[vl->id];
		auto &grow_dir = dk.grow_dir;
		grow_dir.resize(mesh->n_vertices(), -1);
		grow_dir[vl->v.idx()] = 0;
		int nv = mesh->n_vertices();
		int nf = mesh->n_faces();
		std::vector<std::vector<int>> advancing_front[2];
		std::vector<int> hierarchy_vertex[2];
		hierarchy_vertex[0].reserve(pl.size());
		hierarchy_vertex[1].reserve(pl.size());
		BoolVector visited_vl(nv, false);
		visited_vl[vl->v.idx()] = true;
		for (auto& pohl : pl)
		{
			if (!visited_vl[pohl.hl->from / 4])
			{
				hierarchy_vertex[0].push_back(pohl.hl->from / 4);
				grow_dir[pohl.hl->from / 4] = 0;
				visited_vl[pohl.hl->from / 4] = true;
			}
			if (!visited_vl[pohl.hl->to / 4])
			{
				hierarchy_vertex[1].push_back(pohl.hl->to / 4);
				grow_dir[pohl.hl->to / 4] = 1;
				visited_vl[pohl.hl->to / 4] = true;
			}
		}

		for (int i = 0; i < 2; ++i)
		{
			advancing_front[i].push_back(std::move(hierarchy_vertex[i]));
			for (int itertimes = 0; itertimes < 2; ++itertimes)
			{
				const auto& current_af = advancing_front[i].back();
				std::vector<int> hierarchy;
				for (auto caf : current_af)
				{
					for (auto voh : mesh->voh_range(mesh->vertex_handle(caf)))
					{
						VertexHandle vv = voh.to();
						if (visited_vl[vv.idx()])
							continue;
						if (cset.has_face[voh.face().idx()] && cset.has_face[voh.opp().face().idx()])
							continue;
						hierarchy.push_back(vv.idx());
						grow_dir[vv.idx()] = i;
						visited_vl[vv.idx()] = true;
					}
				}
				//if (itertimes == 2)
				{
					for (int id : hierarchy)
					{
						//if (cset.vertex_bound_index[id] != (i == 0 ? dk.from_bound : dk.to_bound))
						//	continue;
						if (cset.vertex_bound_index[id] != dk.from_bound && cset.vertex_bound_index[id] != dk.to_bound)
							continue;
						if (vert_path_similarity_energy[id] > disk_e)
							continue;
						--itertimes;
						break;
					}
				}
				advancing_front[i].push_back(std::move(hierarchy));
			}
		}

		auto &new_vertex = dk.new_vertex;
		auto &newv_flag = dk.new_v_flag;
		auto &region_vertex = dk.region_vertex;
		auto &regionv_flag = dk.region_v_flag;

		newv_flag.resize(nv, false);
		region_vertex.push_back(vl->v);
		regionv_flag.resize(nv, false); regionv_flag[vl->v.idx()] = true;
		for (int i = 0; i < 2; ++i)
		{
			for (auto& ss : advancing_front[i])
			{
				for (auto& tt : ss)
				{
					new_vertex.push_back(mesh->vertex_handle(tt));
					newv_flag[tt] = true;
				}
			}
		}

		auto &new_face = dk.new_face;
		auto &newf_flag = dk.new_f_flag;
		newf_flag.resize(nf, false);
		dk.region_f_flag.resize(nf, false);
		for (auto new_v : new_vertex)
		{
			for (auto vh : mesh->voh_range(new_v))
			{
				if ((newv_flag[vh.to().idx()] || regionv_flag[vh.to().idx()]) &&
					(newv_flag[vh.next().to().idx()] || regionv_flag[vh.next().to().idx()]))
				{
					if (newf_flag[vh.face().idx()] || cset.has_face[vh.face().idx()])
						continue;
					new_face.push_back(vh.face());
					newf_flag[vh.face().idx()] = true;
				}
			}
		}

		//加入边界上x_axis, y_axis的约束
		auto &x_axis = dk.x_axis;
		auto &y_axis = dk.y_axis;
		x_axis.resize(3, mesh->n_faces()); x_axis.setZero();
		y_axis.resize(3, mesh->n_faces()); y_axis.setZero();
		auto &constraint_f_flag = dk.constraint_f_flag;
		constraint_f_flag.resize(mesh->n_faces(), false);
		auto &crossfield = cf->crossfield;
		//可参考《实验记录》2.25示例图
		auto &bound0 = cset.cylinders[dk.from_bound.first].bounds[dk.from_bound.second];
		for (auto hl : bound0)
		{
			int flid = hl->oppo->left;
			if (dk.from_bound.second == 1)
				flid = (flid / 4) * 4 + ((flid + 2) % 4);
			y_axis.col(flid / 4) = crossfield.col(flid);
			x_axis.col(flid / 4) = crossfield.col((flid / 4) * 4 + ((flid + 3) % 4));
			constraint_f_flag[flid / 4] = true;
		}
		auto &bound1 = cset.cylinders[dk.to_bound.first].bounds[dk.to_bound.second];
		for (auto hl : bound1)
		{
			int flid = hl->oppo->left;
			if (dk.to_bound.second == 0)
				flid = (flid / 4) * 4 + ((flid + 2) % 4);
			y_axis.col(flid / 4) = crossfield.col(flid);
			x_axis.col(flid / 4) = crossfield.col((flid / 4) * 4 + ((flid + 3) % 4));
			constraint_f_flag[flid / 4] = true;
		}
	}

	void LoopGen::OptimizeDiskField(disk &dk)
	{
		auto &new_vertex = dk.new_vertex;
		auto &newv_flag = dk.new_v_flag;
		auto &new_face = dk.new_face;
		auto &newf_flag = dk.new_f_flag;
		auto &x_axis = dk.x_axis;
		auto &y_axis = dk.y_axis;

		typedef std::complex<double> COMPLEX;
		std::vector<Eigen::Triplet<COMPLEX>> triple;
		int count = 0;
		std::vector<int> fidmap(mesh->n_faces());
		for (auto f : new_face)
		{
			fidmap[f.idx()] = count++;
		}
		Eigen::VectorXcd b(3 * new_face.size()); b.setZero();
		count = 0;
		auto& position = cf->getPosition();
		auto& faceBase = cf->getFaceBase();
		for (auto f : new_face)
		{
			int fid = f.idx();
			for (auto fh : mesh->fh_range(f))
			{
				if (mesh->is_boundary(fh.edge()))
					continue;
				int gid = fh.opp().face().idx();
				//if()
				if (newf_flag[gid] && fid < gid || dk.constraint_f_flag[gid])
				{
					auto ev = (position.col(fh.to().idx()) - position.col(fh.from().idx())).normalized();
					COMPLEX e_f = COMPLEX(ev.dot(faceBase.col(fid * 2)), -ev.dot(faceBase.col(fid * 2 + 1)));
					COMPLEX e_g = COMPLEX(ev.dot(faceBase.col(gid * 2)), -ev.dot(faceBase.col(gid * 2 + 1)));
					triple.emplace_back(count, fidmap[fid], e_f);
					if (newf_flag[gid])
					{
						triple.emplace_back(count, fidmap[gid], -e_g);
					}
					else
					{
						COMPLEX dir = COMPLEX(x_axis.col(gid).dot(faceBase.col(2 * gid)), x_axis.col(gid).dot(faceBase.col(2 * gid + 1)));
						b(count) += e_g * dir;
					}
					++count;
				}
			}
		}

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> solver;
		Eigen::SparseMatrix<COMPLEX> A(count, new_face.size());
		A.setFromTriplets(triple.begin(), triple.end());
		Eigen::SparseMatrix<COMPLEX> AH = A.adjoint();
		solver.compute(AH*A);
		b = AH * b.head(count);
		b = solver.solve(b);
		for (auto f : new_face)
		{
			int fid = f.idx();
			double theta = std::arg(b(fidmap[fid]));
			x_axis.col(fid) = faceBase.col(fid * 2) * cos(theta) + faceBase.col(fid * 2 + 1) * sin(theta);
			y_axis.col(fid) = faceBase.col(fid * 2) * cos(theta + 0.5 * PI) + faceBase.col(fid * 2 + 1) * sin(theta + 0.5 * PI);
		}
	}

	bool LoopGen::RefinePathByField(VertexHandle v, disk &dk, BoolVector &visited_v, BoolVector &visited_f)
	{
		Eigen::Vector4d plane; plane.setZero();
		double fromdis = 0; double todis = 0;
		auto sign_dist = [&](OpenMesh::Vec3d &p)
		{
			return plane(0)*p[0] + plane(1)*p[1] + plane(2)*p[2] + plane(3);
		};
		PlaneLoop field_path;

		auto &y_axis = dk.y_axis;
		if (cset.vertex_bound_index[v.idx()] != dk.from_bound)
		{
			HalfedgeHandle hh; hh.invalidate();
			plane.setZero();
			for (auto vf : mesh->vf_range(v))
				plane.head(3) += y_axis.col(vf.idx());
			plane.head(3).normalize();
			for (auto voh : mesh->voh_range(v))
			{
				if (!visited_f[voh.face().idx()])
					continue;
				//plane.head(3) = y_axis.col(voh.face().idx());
				auto &pos = mesh->point(v);
				plane(3) = -(pos[0] * plane(0) + pos[1] * plane(1) + pos[2] * plane(2));
				fromdis = sign_dist(mesh->point(voh.next().to()));
				todis = sign_dist(mesh->point(voh.to()));
				if (fromdis <= 0 && todis >= 0)
				{
					hh = voh.next().opp();
					break;
				}
			}
			if (!hh.is_valid())
				return false;
			double c = todis / (todis - fromdis);
			field_path.emplace_back(&m4.halfedgelayers[hh.idx() * 4], std::min(0.99, std::max(0.01, c)));
			//while (!cset.bound_edge_flag[hh.idx() / 2])
			while (cset.vertex_bound_index[mesh->to_vertex_handle(hh).idx()] != dk.from_bound || !cset.bound_edge_flag[hh.idx() / 2])
			{
				if (!visited_f[mesh->face_handle(hh).idx()])
					return false;
				plane.head(3) = y_axis.col(mesh->face_handle(hh).idx());
				auto pos = field_path.back().point(m4);
				plane(3) = -(pos[0] * plane(0) + pos[1] * plane(1) + pos[2] * plane(2));
				double dis = sign_dist(mesh->point(mesh->to_vertex_handle(mesh->next_halfedge_handle(hh))));
				if (dis > 0)
				{
					hh = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(hh));
					todis = dis;
					fromdis = sign_dist(mesh->point(mesh->from_vertex_handle(hh)));
					c = todis / (todis - fromdis);
					field_path.emplace_back(&m4.halfedgelayers[hh.idx() * 4], std::min(0.99, std::max(0.01, c)));
				}
				else
				{
					hh = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(hh));
					fromdis = dis;
					todis = sign_dist(mesh->point(mesh->to_vertex_handle(hh)));
					c = todis / (todis - fromdis);
					field_path.emplace_back(&m4.halfedgelayers[hh.idx() * 4], std::min(0.99, std::max(0.01, c)));
				}
			}
		}
		std::reverse(field_path.begin(), field_path.end());
		field_path.emplace_back(nullptr, v.idx()*1.0);
		if (cset.vertex_bound_index[v.idx()] != dk.to_bound)
		{
			HalfedgeHandle hh; hh.invalidate();
			plane.setZero();
			for (auto vf : mesh->vf_range(v))
				plane.head(3) += y_axis.col(vf.idx());
			plane.head(3).normalize();
			for (auto voh : mesh->voh_range(v))
			{
				if (!visited_f[voh.face().idx()])
					continue;
				//plane.head(3) = y_axis.col(voh.face().idx());
				auto &pos = mesh->point(v);
				plane(3) = -(pos[0] * plane(0) + pos[1] * plane(1) + pos[2] * plane(2));
				fromdis = sign_dist(mesh->point(voh.to()));
				todis = sign_dist(mesh->point(voh.next().to()));
				if (fromdis <= 0 && todis >= 0)
				{
					hh = voh.next();
					break;
				}
			}
			if (!hh.is_valid())
				return false;
			double c = todis / (todis - fromdis);
			field_path.emplace_back(&m4.halfedgelayers[hh.idx() * 4], std::min(0.99, std::max(0.01, c)));
			while (cset.vertex_bound_index[mesh->to_vertex_handle(hh).idx()] != dk.to_bound || !cset.bound_edge_flag[hh.idx() / 2])
			{
				hh = mesh->opposite_halfedge_handle(hh);
				if (!visited_f[mesh->face_handle(hh).idx()])
					return false;
				plane.head(3) = y_axis.col(mesh->face_handle(hh).idx());
				auto pos = field_path.back().point(m4);
				plane(3) = -(pos[0] * plane(0) + pos[1] * plane(1) + pos[2] * plane(2));
				double dis = sign_dist(mesh->point(mesh->to_vertex_handle(mesh->next_halfedge_handle(hh))));
				if (dis > 0)
				{
					hh = mesh->next_halfedge_handle(hh);
					todis = dis;
					fromdis = sign_dist(mesh->point(mesh->from_vertex_handle(hh)));
					c = todis / (todis - fromdis);
					field_path.emplace_back(&m4.halfedgelayers[hh.idx() * 4], std::min(0.99, std::max(0.01, c)));
				}
				else
				{
					hh = mesh->prev_halfedge_handle(hh);
					fromdis = dis;
					todis = sign_dist(mesh->point(mesh->to_vertex_handle(hh)));
					c = todis / (todis - fromdis);
					field_path.emplace_back(&m4.halfedgelayers[hh.idx() * 4], std::min(0.99, std::max(0.01, c)));
				}
			}
		}
		all_vertice_path[v.idx()] = std::move(field_path);
		return true;
	}

	double LoopGen::EvaluatePathSimilarity(Eigen::Matrix3Xd &path0, Eigen::Matrix3Xd &path1)
	{
		//对两个path重新采样，对比采样点的切向，从而定义相似性
		int n = path0.cols() + path1.cols();
		/*static int sefz = 0;
		++sefz;*/
		auto pathLength = [&](Eigen::Matrix3Xd& path, Eigen::VectorXd& seg)
		{
			int cols = path.cols() - 1;
			seg.resize(cols); seg.setZero();
			double sum = 0;
			for (int i = 0; i < cols; ++i)
			{
				seg(i) = (path.col(i + 1) - path.col(i)).norm();
				sum += seg(i);
			}
			return sum;
		};
		auto assembleFragment = [&](Eigen::Matrix3Xd& fragment, Eigen::Matrix3Xd& path)
		{
			int cols = path.cols();
			Eigen::VectorXd seg;
			double step = pathLength(path, seg) / n;
			Eigen::Vector3d vec = (path.col(1) - path.col(0)).normalized();
			fragment.col(0) = vec;
			int r = 0;
			double l = seg(0);
			int mark = 0;
			double u0 = 0;
			for (int i = 1; i < n; ++i)
			{
				u0 += step;
				if (u0 > l)
				{
					while (u0 > l)
					{
						++r; ++mark;
						l += seg(r % cols);
					}
					vec = (path.col((mark + 1) % cols) - path.col(mark % cols)).normalized();
					//if(sefz==438)
					//dprint("mark", mark, cols);
					/*if (sefz == 438)
					{
						dprint(i, u0, l, mark, r);
						dprint(path(0, (mark + 1) % cols), path(1, (mark + 1) % cols), path(2, (mark + 1) % cols));
						dprint(path(0, mark%cols), path(1, mark%cols), path(2, mark%cols));
					}*/
				}
				fragment.col(i) = vec;
				/*if (sefz == 438)
				{
					dprint("faeifhaefjs", i, vec(0),vec(1),vec(2));
				}*/
			}
		};
		/*if (sefz == 438)
		{
			for (int i = 0; i < path0.cols(); ++i)
				dprint(i, path0(0, i), path0(1, i), path0(2, i));
			for (int i = 0; i < path1.cols(); ++i)
				dprint(i, path1(0, i), path1(1, i), path1(2, i));
		}*/
		if (path0.cols() < 2 || path1.cols() < 2)
			return YYSS_INFINITE;
		Eigen::Matrix3Xd fragment0, fragment1;
		fragment0.resize(3, n); fragment1.resize(3, n);
		assembleFragment(fragment0, path0);
		assembleFragment(fragment1, path1);
		double sum = 0;
		double dot = 0;
		//dprint(sefz);
		for (int i = 0; i < n; ++i)
		{
			/*if (sefz == 438)
			{
				dprint(i);
				dprint(fragment0(0, i), fragment0(1, i), fragment0(2, i));
				dprint(fragment1(0, i), fragment1(1, i), fragment1(2, i));
			}*/
			for (int j = i + 1; j < n; ++j)
			{
				double theta = conservativeArcCos(fragment0.col(i).dot(fragment0.col(j))) - conservativeArcCos(fragment1.col(i).dot(fragment1.col(j)));
				sum += fabs(theta);
				//if (sefz == 438)
				//{
				//	//dprint(i,j,theta);
				//	//dprint(fragment0(0,i),fragment0(1,i),fragment0(2,i),)
				//}
			}
		}
		/*if (sefz == 438)
		{
			int p = 0;
		}*/
		return 2.0 * sum / (n * (n - 1));
	}

	void LoopGen::AssembleSampling(PlaneLoop &path, Eigen::Matrix3Xd &sampling, int path_sampling_num)
	{
		Eigen::Matrix3Xd path_pos;
		path_pos.resize(3, path.size()); path_pos.setZero();
		for (int i = 0; i < path.size(); ++i)
		{
			auto pos = path[i].point(m4);
			path_pos.col(i) << pos[0], pos[1], pos[2];
		}
		Eigen::Vector3d barycenter;
		auto pathLength = [&](Eigen::VectorXd& seg)
		{
			int cols = path_pos.cols() - 1;
			seg.resize(cols); seg.setZero();
			barycenter.setZero();
			double sum = 0;
			for (int i = 0; i < cols; ++i)
			{
				seg(i) = (path_pos.col(i + 1) - path_pos.col(i)).norm();
				sum += seg(i);
				barycenter += path_pos.col(i);
			}
			barycenter = (barycenter + path_pos.col(cols)) / (cols + 1);
			return sum;
		};

		int cols = path.size();
		Eigen::VectorXd seg;
		double len = pathLength(seg);
		double step = len / path_sampling_num;
		sampling.resize(3, path_sampling_num); sampling.setZero();
		sampling.col(0) = (path_pos.col(0) - barycenter) / len;
		int r = 0;
		double l = seg(0);
		int mark = 0;
		double u0 = 0;
		Eigen::Vector3d vec = (path_pos.col(1) - path_pos.col(0)).normalized();
		for (int i = 1; i < path_sampling_num - 1; ++i)
		{
			u0 += step;
			if (u0 > l)
			{
				while (u0 > l)
				{
					++r; ++mark;
					l += seg(r % cols);
				}
				vec = (path_pos.col((mark + 1) % cols) - path_pos.col(mark % cols)).normalized();
			}
			//sampling.col(i)=(path_pos.col(mark%cols)+u0*)//github的第66行，这个u0是不是有问题？为什么是u0*vec，感觉不太对
			sampling.col(i) = (path_pos.col((mark + 1) % cols) + (u0 - l)*vec - barycenter) / len;
		}
		sampling.col(path_sampling_num - 1) = (path_pos.col(cols - 1) - barycenter) / len;
	}

	double LoopGen::AssembleSimilarityAngle(PlaneLoop &path, Eigen::VectorXd &sa, int path_fragment_num)
	{
		Eigen::Matrix3Xd path_pos;
		path_pos.resize(3, path.size()); path_pos.setZero();
		for (int i = 0; i < path.size(); ++i)
		{
			auto pos = path[i].point(m4);
			path_pos.col(i) << pos[0], pos[1], pos[2];
		}
		auto pathLength = [&](Eigen::VectorXd& seg)
		{
			int cols = path_pos.cols() - 1;
			seg.resize(cols); seg.setZero();
			double sum = 0;
			for (int i = 0; i < cols; ++i)
			{
				seg(i) = (path_pos.col(i + 1) - path_pos.col(i)).norm();
				sum += seg(i);
			}
			return sum;
		};

		Eigen::Matrix3Xd fragment;
		fragment.resize(3, path_fragment_num); fragment.setZero();
		int cols = path.size();
		Eigen::VectorXd seg;
		double step = pathLength(seg) / path_fragment_num;
		Eigen::Vector3d vec = (path_pos.col(1) - path_pos.col(0)).normalized();
		fragment.col(0) = vec;
		int r = 0;
		double l = seg(0);
		int mark = 0;
		double u0 = 0;
		for (int i = 1; i < path_fragment_num; ++i)
		{
			u0 += step;
			if (u0 > l)
			{
				while (u0 > l)
				{
					++r; ++mark;
					l += seg(r % cols);
				}
				vec = (path_pos.col((mark + 1) % cols) - path_pos.col(mark % cols)).normalized();
			}
			fragment.col(i) = vec;
		}
		sa.resize(path_fragment_num*(path_fragment_num + 1) / 2); sa.setZero();
		int count = 0;
		for (int i = 0; i < path_fragment_num; ++i)
			for (int j = i + 1; j < path_fragment_num; ++j)
			{
				sa(count++) = conservativeArcCos(fragment.col(i).dot(fragment.col(j)));
				/*if (ie == 124)
					dprint(i, j, sa(count - 1));*/
			}
		/*if (ie == 124)
		{
			int lp = 0;
		}*/
		return step * path_fragment_num;
	}

	bool LoopGen::SpreadSubRegion(disk &dk, bool grow_flag[2])
	{
		int nv = mesh->n_vertices();
		int nf = mesh->n_faces();

		auto& new_vertex = dk.new_vertex;
		auto& new_face = dk.new_face;
		auto &region_vertex = dk.region_vertex;
		auto &region_face = dk.region_face;
		auto& newv_flag = dk.new_v_flag;
		auto& newf_flag = dk.new_f_flag;
		auto &regionv_flag = dk.region_v_flag;
		auto &regionf_flag = dk.region_f_flag;
		BoolVector visited_v = dk.region_v_flag;
		BoolVector visited_f = dk.region_f_flag;
		for (auto new_v : new_vertex)
		{
			visited_v[new_v.idx()] = true;
			newv_flag[new_v.idx()] = false;
		}
		for (auto new_f : new_face)
		{
			visited_f[new_f.idx()] = true;
			newf_flag[new_f.idx()] = false;
		}
		//dprint("new vertex:", lp.new_vertex.size());
		//提取含有完整loop的区域

		if (region_vertex.size() == 1)
			if (!RefinePathByField(region_vertex.front(), dk, visited_v, visited_f))
				return false;
		std::vector<VertexHandle> vertex_cache;
		BoolVector vertex_cache_flag = dk.region_v_flag;
		//return true;
		for (auto new_v : new_vertex)
		{
			//dprint(new_v.idx());
			if (new_v.idx() == 13397 && region_vertex.front().idx() == 23937/4)
			{
				int p = 0;
			}
			if (RefinePathByField(new_v, dk, visited_v, visited_f))
			{
				vertex_cache_flag[new_v.idx()] = true;
				vertex_cache.push_back(new_v);
			}
		}
		//dprint("vertex cache:", vertex_cache.size());
		v_cache_flag = vertex_cache_flag;
#if PRINT_DEBUG_INFO
		dprint("提取含有完整loop的区域");
#endif
		/*if (region_face.empty())
		{
			auto hl_begin = region_vertex.front()->hl;
			auto hl_transfer = hl_begin;
			do
			{
				if (!vertex_cache_flag[hl_transfer->to])
					return false;
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			if (!CheckTopology(vertex_cache, vertex_cache_flag, lp.grow_dir))
				return false;
		}*/
		//若区域不连接两个边界，则退出
		bool fromflag = false; bool toflag = false;
		for (auto v : vertex_cache)
		{
			if (cset.vertex_bound_index[v.idx()] == dk.from_bound)
				fromflag = true;
			if (cset.vertex_bound_index[v.idx()] == dk.to_bound)
				toflag = true;
		}
		if (!fromflag || !toflag)
			return false;
		//若区域不满足disk拓扑，则退出
		if (region_face.empty())
			if (!CheckDiskTopology(vertex_cache, vertex_cache_flag))
				return false;
		if (vertex_cache.empty())
			return false;

#if PRINT_DEBUG_INFO
		dprint("拓扑检查");
#endif

		if (region_face.empty())
		{
			Eigen::Vector3d from_dir; from_dir.setZero();
			Eigen::Vector3d to_dir; to_dir.setZero();
			//for (auto v : vertex_cache)
			int min_id = -1;
			double min_e = YYSS_INFINITE;
			for (int k = 0; k < vertex_cache.size(); ++k)
			{
				auto &path = all_vertice_path[vertex_cache[k].idx()];
				from_dir.setZero();
				to_dir.setZero();
				if (path.front().hl)
				{
					from_dir = dk.y_axis.col(path.front().hl->left / 4) + dk.y_axis.col(path.front().hl->oppo->left / 4);
					from_dir.normalize();
				}
				else
				{
					for (auto vf : mesh->vf_range(mesh->vertex_handle(path.front().c)))
						from_dir += dk.y_axis.col(vf.idx());
					from_dir.normalize();
				}
				if (path.back().hl)
				{
					to_dir = dk.y_axis.col(path.back().hl->left / 4) + dk.y_axis.col(path.back().hl->oppo->left / 4);
					to_dir.normalize();
				}
				else
				{
					for (auto vf : mesh->vf_range(mesh->vertex_handle(path.back().c)))
						to_dir += dk.y_axis.col(vf.idx());
					to_dir.normalize();
				}
				double e = 0;
				double length = 0;
				Eigen::Vector3d pos;
				OpenMesh::Vec3d knot[2];
				knot[0] = path.front().point(m4);
				for (int i = 1; i < path.size(); ++i)
				{
					knot[i % 2] = path[i].point(m4);
					for (int j = 0; j < 3; ++j)
						pos(j) = knot[i % 2][j] - knot[(i + 1) % 2][j];
					e += fabs(pos.dot(from_dir)) + fabs(pos.dot(to_dir));
					length += pos.norm();
				}
				e /= length * 2.0;
				if (e < min_e)
				{
					min_id = k;
					min_e = e;
				}
			}

			if (min_id != -1)
			{
#if 1
				auto v = region_vertex.front();
				region_vertex.front() = vertex_cache[min_id];
				dk.grow_dir[v.idx()] = 0;
				regionv_flag[vertex_cache[min_id].idx()] = true;
				vertex_cache[min_id] = v;
				regionv_flag[v.idx()] = false;
#endif
			}
		}

		auto& grow_dir = dk.grow_dir;

		//检测相似性能量
#if USE_NEW_SIMILARITY_ENERGY
		auto &normal_sampling = dk.normal_sampling;
		int path_sampling_num = all_vertice_path[region_vertex.front().idx()].size();
		if (path_sampling_num == 0)
			return false;
		if (vertex_cache_flag[region_vertex.front().idx()] && !dk.has_ns)
		{
			dk.has_ns = true;
			AssembleSampling(all_vertice_path[region_vertex.front().idx()], normal_sampling, path_sampling_num);
		}
#else
		auto& normal_similarity_angle = dk.normal_similarity_angle;
		int path_fragment_num = all_vertice_path[region_vertex.front().idx()].size();
		if (path_fragment_num == 0)
			return false;
		if (vertex_cache_flag[region_vertex.front().idx()] && !dk.has_nsa)
		{
			dk.has_nsa = true;
			//auto pohl = all_vertice_path[region_vertex.front().idx()].back();
			//double len = AssembleSimilarityAngle(region_vertex.front(), normal_similarity_angle, tn, path_fragment_num);
			double len = AssembleSimilarityAngle(all_vertice_path[region_vertex.front().idx()], normal_similarity_angle, path_fragment_num);
			dk.length[0] = std::min(dk.length[0], len);
			dk.length[1] = std::max(dk.length[1], len);
			dk.length[2] = std::max(dk.length[2], len);
		}
#endif

		BoolVector if_similarity_energy_low;
		//if (region_vertex.size() > 1)
		{
			if_similarity_energy_low.resize(nv, false);
			int exceed[2] = { 0,0 };
			//#pragma omp parallel for
			for (int i = 0; i < vertex_cache.size(); ++i)
			{
				int new_id = vertex_cache[i].idx();
				int grow_id = grow_dir[new_id];
				if (!grow_flag[grow_id])
					continue;
#if USE_NEW_SIMILARITY_ENERGY
				Eigen::Matrix3Xd sampling;
				AssembleSampling(all_vertice_path[vertex_cache[i].idx()], sampling, path_sampling_num);
				double ey = ICP_Energy(normal_sampling, sampling);
				if (ey < disk_e)
					if_similarity_energy_low[new_id] = true;
#else
				Eigen::VectorXd similarity_angle;
				//double len = AssembleSimilarityAngle(vertex_cache[i], similarity_angle, tn, path_fragment_num);
				double len = AssembleSimilarityAngle(all_vertice_path[vertex_cache[i].idx()], similarity_angle, path_fragment_num);
				double sum = 0;
				//double seaa = 0;
				for (int j = 0; j < similarity_angle.size(); ++j) {
					sum += fabs(normal_similarity_angle(j) - similarity_angle(j));
					//seaa = std::max(seaa, fabs(normal_similarity_angle(j) - similarity_angle(j)));
					if (region_vertex.front().idx() == 147012 / 4 && i == 463)
					{
						//dprint(i, j, normal_similarity_angle(j), similarity_angle(j), fabs(normal_similarity_angle(j) - similarity_angle(j)));
					}
				}
				dk.length[0] = std::min(dk.length[0], len);
				dk.length[grow_id + 1] = std::max(dk.length[grow_id + 1], len);
				//dprint(i, sum / similarity_angle.size());
				if (sum < disk_e * similarity_angle.size()/* && tn.length[0] * 1.3 > tn.length[grow_id + 1]*/)
				{
					if_similarity_energy_low[new_id] = true;
				}
#endif
			}
			/*if (region_vertex.front().idx() == 147012 / 4)
			{
				grow_flag[0] = false; grow_flag[1] = false;
				newf_flag = visited_f;
				return false;
			}*/
			//if_similarity_energy_low.resize(nv, true);
#if PRINT_DEBUG_INFO
			dprint("检测相似性能量");
#endif
			for (auto vc : vertex_cache)
			{
				if (!grow_flag[grow_dir[vc.idx()]])
				{
					if_similarity_energy_low[vc.idx()] = false;
					continue;
				}
				if (!if_similarity_energy_low[vc.idx()])
					++exceed[grow_dir[vc.idx()]];
				if (exceed[grow_dir[vc.idx()]] > 1)
					grow_flag[grow_dir[vc.idx()]] = false;
			}
		}
		//else
		//	if_similarity_energy_low.resize(nv, true);
		//else
			//if_similarity_energy_low.resize(nv, true);


		//更新新区域
		int count = region_vertex.size();
		int begin_ = count;
		for (auto new_v : vertex_cache)
		{
			if (if_similarity_energy_low[new_v.idx()])
			{
				region_vertex.push_back(new_v);
				int newid = new_v.idx();
				regionv_flag[newid] = true;
			}
		}

		for (auto new_f : new_face)
		{
			for (auto fv : mesh->fv_range(new_f))
			{
				if (!regionv_flag[fv.idx()])
					goto target1;
			}
			regionf_flag[new_f.idx()] = true;
			region_face.push_back(new_f);
		target1:;
		}

		//若区域不连接两个边界，则退出
		fromflag = false; toflag = false;
		for (auto v : region_vertex)
		{
			if (cset.vertex_bound_index[v.idx()] == dk.from_bound)
				fromflag = true;
			if (cset.vertex_bound_index[v.idx()] == dk.to_bound)
				toflag = true;
		}
		if (!fromflag || !toflag)
		{
			region_vertex.clear();
			region_face.clear();
			regionv_flag.resize(nv, false);
			regionf_flag.resize(nf, false);
			return false;
		}
		//若区域不满足disk拓扑，则退出
		if (region_face.empty())
			return false;
		if (!CheckDiskTopology(region_vertex, regionv_flag))
		{
			region_vertex.clear();
			region_face.clear();
			regionv_flag.resize(nv, false);
			regionf_flag.resize(nf, false);
			return false;
		}

#if 1
#if PRINT_DEBUG_INFO
		dprint("更新新区域");
		dprint("当前区域含有的顶点数：", region_vertex.size());
#endif
		//若已超过能量阈值或者遇到接近平面的区域，则退出
		if (!grow_flag[0] && !grow_flag[1])
			return false;

		//扩展advancing_front
		newv_flag.resize(nv, false);
		new_vertex.clear();
		for (int i = begin_; i < count; ++i)
		{
			auto rvi = region_vertex[i];
			int growid = grow_dir[rvi.idx()];
			if (!grow_flag[growid])
				continue;
			/*for (auto vv : mesh->vv_range(rvi))
			{
				int vvid = vv.idx();
				if (regionv_flag[vvid] || newv_flag[vvid])
					continue;
				if (cset.has_vertex[vvid] && cset.vertex_bound_index[vvid].first != tn.from_bound.first
					&& cset.vertex_bound_index[vvid].first != tn.to_bound.first)
					continue;
				new_vertex.push_back(vv);
				newv_flag[vvid] = true;
				grow_dir[vvid] = growid;
			}*/
			for (auto voh : mesh->voh_range(rvi))
			{
				VertexHandle vv = voh.to();
				if (regionv_flag[vv.idx()] || newv_flag[vv.idx()])
					continue;
				if (cset.has_face[voh.face().idx()] && cset.has_face[voh.opp().face().idx()])
					continue;
				new_vertex.push_back(vv);
				newv_flag[vv.idx()] = true;
				grow_dir[vv.idx()] = growid;
			}
		}
		begin_ = 0;
		for (int i = 0; i < extend_layer - 1; ++i)
		{
			int end_ = new_vertex.size();
			for (int j = begin_; j < end_; ++j)
			{
				auto rvj = new_vertex[j];
				int growid = grow_dir[rvj.idx()];
				/*for (auto vv : mesh->vv_range(rvj))
				{
					int vvid = vv.idx();
					if (regionv_flag[vvid] || newv_flag[vvid])
						continue;
					if (cset.has_vertex[vvid] && cset.vertex_bound_index[vvid].first != tn.from_bound.first
						&& cset.vertex_bound_index[vvid].first != tn.to_bound.first)
						continue;
					new_vertex.push_back(vv);
					newv_flag[vvid] = true;
					grow_dir[vvid] = growid;
				}*/
				for (auto voh : mesh->voh_range(rvj))
				{
					VertexHandle vv = voh.to();
					if (regionv_flag[vv.idx()] || newv_flag[vv.idx()])
						continue;
					if (cset.has_face[voh.face().idx()] && cset.has_face[voh.opp().face().idx()])
						continue;
					new_vertex.push_back(vv);
					newv_flag[vv.idx()] = true;
					grow_dir[vv.idx()] = growid;
				}
			}
			begin_ = end_;
		}
		if (new_vertex.empty())
			return false;
		newf_flag.resize(nf, false);
		new_face.clear();
		for (auto new_v : new_vertex)
		{
			int newid = new_v.idx();
			for (auto voh : mesh->voh_range(new_v))
			{
				int vfid = voh.face().idx();
				if (regionf_flag[vfid] || newf_flag[vfid] || cset.has_face[vfid])
					continue;
				if (!regionv_flag[voh.to().idx()] && !newv_flag[voh.to().idx()])
					continue;
				if (!regionv_flag[voh.next().to().idx()] && !newv_flag[voh.next().to().idx()])
					continue;
				new_face.push_back(voh.face());
				newf_flag[voh.face().idx()] = true;
			}
		}
		if (new_face.empty())
			return false;
#if PRINT_DEBUG_INFO
		dprint("扩展advancing_front");
#endif
		//if (region_vertex.front().idx() == 147012 / 4)
		//	return false;
		//优化新区域的场
		for (auto rf : region_face)
			dk.constraint_f_flag[rf.idx()] = true;
		//set_one_field(tn);
		//OptimizeDiskField(dk);
#if PRINT_DEBUG_INFO
		dprint("优化新区域的场");
#endif
#endif
		return true;
	}

	bool LoopGen::CheckDiskTopology(std::vector<VertexHandle> &vertex_set, BoolVector &vs_flag)
	{
		if (vertex_set.empty())
			return false;

#pragma region initialize fs_flag
		int nv = mesh->n_vertices();
		int nf = mesh->n_faces();
		BoolVector fs_flag(nf, false);
		for (auto v : vertex_set)
		{
			for (auto voh : mesh->voh_range(v))
			{
				if (cset.has_face[voh.face().idx()])
					continue;
				if (!vs_flag[voh.to().idx()] || !vs_flag[voh.next().to().idx()])
					continue;
				fs_flag[voh.face().idx()] = true;
			}
		}
#pragma endregion

#pragma region check connectivity
		BoolVector visited(nv, false);
		std::queue<VertexHandle> tree;
		tree.push(vertex_set.front());
		visited[vertex_set.front().idx()] = true;
		int count = 0;
		while (!tree.empty())
		{
			auto v = tree.front();
			tree.pop();
			++count;
			for (auto vv : mesh->vv_range(v))
			{
				if (visited[vv.idx()] || !vs_flag[vv.idx()])
					continue;
				auto hh = mesh->find_halfedge(v, vv);
				if (!fs_flag[hh.face().idx()] && !fs_flag[hh.opp().face().idx()])
					continue;
				tree.push(vv);
				visited[vv.idx()] = true;
			}
		}
		if (count < vertex_set.size())
			return false;
#pragma endregion

#pragma region check manifold
		BoolVector bv_flag(nv, false);
		for (auto v : vertex_set)
		{
			for (auto vf : mesh->vf_range(v))
			{
				if (!fs_flag[vf.idx()])
				{
					bv_flag[v.idx()] = true;
					break;
				}
			}
		}
#pragma endregion

#pragma region check boundary number
		{
			BoolVector be_flag(mesh->n_edges(), false);
			/*static int es = 0;
			++es;
			dprint(es);*/
			count = 0;
			for (auto v : vertex_set)
			{
				/*if (es == 40 && v.idx() == 1942)
				{
					dprint(v.idx());
				}*/
				if (!bv_flag[v.idx()])
					continue;
				/*auto he = mesh->voh_begin(v).handle();
				while (fs_flag[mesh->face_handle(he).idx()] || !fs_flag[mesh->opposite_face_handle(he).idx()])
					he = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));*/
				HalfedgeHandle he; he.invalidate();
				for (auto voh : mesh->voh_range(v))
				{
					if (!fs_flag[mesh->face_handle(voh).idx()] && fs_flag[mesh->opposite_face_handle(voh).idx()])
					{
						he = voh; break;
					}
				}
				if (!he.is_valid())
					return false;
				if (be_flag[he.idx() / 2])
					continue;
				HalfedgeHandle he_transfer = he;
				do
				{
					be_flag[he_transfer.idx() / 2] = true;
					he_transfer = mesh->opposite_halfedge_handle(he_transfer);
					while (fs_flag[mesh->face_handle(he_transfer).idx()] || !fs_flag[mesh->opposite_face_handle(he_transfer).idx()])
						he_transfer = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(he_transfer));
				} while (he_transfer != he);
				++count;
				if (count > 1)
					return false;
			}
			if (count < 1)
				return false;
		}
#pragma endregion
		return true;
	}

}