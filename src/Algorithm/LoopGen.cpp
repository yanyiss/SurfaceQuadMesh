#include "LoopGen.h"
#include <omp.h>

#define COMPUTE_NEW_PLANELOOP 1
#define COMPUTE_NEW_ENERGY 1

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

	double LoopGen::EvaluateSimilarity(Eigen::Matrix3Xd& loop0, Eigen::Matrix3Xd& loop1, double u, int begin_seg)
	{
		//对两个loop重新采样，对比采样点的切向，从而定义相似性
		int n = loop0.cols() + loop1.cols();
		auto loopLength = [&](Eigen::Matrix3Xd& loop, Eigen::VectorXd& seg, int mark)
		{
			int cols = loop.cols();
			double sum = 0;
			for (int i = 0; i < cols; ++i)
			{
				seg(i) = (loop.col((mark + i + 1) % cols) - loop.col((mark + i) % cols)).norm();
				sum += seg(i);
			}
			return sum;
		};
		auto assembleFragment = [&](Eigen::Matrix3Xd& fragment, double u0, int mark, Eigen::Matrix3Xd& loop)
		{
			int cols = loop.cols();
			Eigen::VectorXd seg(cols);
			double step = loopLength(loop, seg, mark) / n;
			Eigen::Vector3d vec = (loop.col((mark + 1) % cols) - loop.col(mark % cols)).normalized();
			fragment.col(0) = vec;
			int r = 0;
			double l = seg(0);
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
				fragment.col(i) = vec;
			}
		};

		Eigen::Matrix3Xd fragment0, fragment1;
		fragment0.resize(3, n); fragment1.resize(3, n);
		assembleFragment(fragment0, 0, 0, loop0);
		assembleFragment(fragment1, u, begin_seg, loop1);
		double sum = 0;
		double dot = 0;
		for (int i = 0; i < n; ++i)
		{
			for (int j = i + 1; j < n; ++j)
			{
				double theta = conservativeArcCos(fragment0.col(i).dot(fragment0.col(j))) - conservativeArcCos(fragment1.col(i).dot(fragment1.col(j)));
				sum += fabs(theta);
			}
		}
		return 2.0 * sum / (n * (n - 1));
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

	bool LoopGen::FieldAligned_PlanarLoop(/*M4 &m4, */VertexLayer* vl, std::vector<VertexLayer*> &loop)
	{
		struct LayerNode
		{
			int id;
			int count;
			double dist;
			LayerNode() {}
			LayerNode(int id_, double dist_, int count_) :id(id_), dist(dist_), count(count_) {}
			bool operator>(const LayerNode& x) const { return dist > x.dist; }
		};
		std::priority_queue<LayerNode, std::vector<LayerNode>, std::greater<LayerNode>> pq;

		int nvl = m4.verticelayers.size();
		int vid;
		std::vector<double> distance(nvl, YYSS_INFINITE);
		std::vector<int> count(nvl, 0);
		std::vector<HalfedgeLayer*> prev(nvl, nullptr);
		BoolVector visited(nvl, false);

		HalfedgeLayer* hl_begin = vl->hl;
		HalfedgeLayer* hl_transfer = hl_begin;
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
			if (fromid == vl->id)
			{
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

		loop.clear();
		loop.push_back(vl);
		VertexLayer* prev_layer = &m4.verticelayers[prev[vl->id]->from];
		while (prev_layer != vl)
		{
			loop.push_back(prev_layer);
			prev_layer = &m4.verticelayers[prev[prev_layer->id]->from];
		}
		loop.push_back(vl);
		std::reverse(loop.begin(), loop.end());
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

	void LoopGen::InitializePQ()
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
			/*if (m4.sing_flag[m4.verticelayers[i].v.idx()])
			{
				InfoOnMesh[i].dir = Forward;
				InfoOnMesh[i].plid = count;
				++count;
			}
			else
			{
				int type = i - m4.verticemap[m4.verticelayers[i].v.idx()];
				switch (type)
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
			}*/
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
				if (FieldAligned_PlanarLoop(vl, loop))
					RefineLoopByPlanarity(loop, pls[InfoOnMesh[vl->id].plid]);
				++vl;
				if (FieldAligned_PlanarLoop(vl, loop))
					RefineLoopByPlanarity(loop, pls[InfoOnMesh[vl->id].plid]);
			}
			WritePlaneLoop(pls, model_name, mesh);
		}
		tr.out("Time of Initializing Planar Loops on All Vertices:");

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
					similarity_energy[hlid] = EvaluateSimilarity(loop0, loop1, u0, id);
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
						if (InfoOnMesh[toid].energy > energy_threshold)
							goto target0;
						visited_vl[toid] = true;
						hierarchy.push_back(toid);
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);
			}
			advancing_front[0].push_back(std::move(hierarchy));
		}
	target0:;
		while (true)
		{
			const auto& current_af = advancing_front[1].back();
			std::vector<int> hierarchy;
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
						if (InfoOnMesh[toid].energy > energy_threshold)
							goto target1;
						visited_vl[toid] = true;
						hierarchy.push_back(toid);
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);
			}
			advancing_front[1].push_back(std::move(hierarchy));
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
					return false;
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			if (!CheckTopology(vertex_cache, vertex_cache_flag, lp.grow_dir))
				return false;
		}
		if (vertex_cache.empty())
			return false;

#if PRINT_DEBUG_INFO
		dprint("拓扑检查");
#endif

		auto& regionv_flag = lp.region_v_flag;
		auto& regionf_flag = lp.region_f_flag;
		auto& grow_dir = lp.grow_dir;

		//检测相似性能量
		auto& normal_similarity_angle = lp.normal_similarity_angle;
		auto& all_pl = lp.all_pl;
		int loop_fragment_num = all_pl[lp.region_vertex.front()->id].size();
		if (vertex_cache_flag[lp.region_vertex.front()->id] && !lp.has_nsa)
		{
			lp.has_nsa = true;
			AssembleSimilarityAngle(lp.region_vertex.front(), normal_similarity_angle, lp, loop_fragment_num);
		}

		BoolVector if_similarity_energy_low;
		if (lp.region_vertex.size() > 1)
		{
			if_similarity_energy_low.resize(nvl, false);
			int exceed[2] = { 0,0 };
#pragma omp parallel for
			for (int i = 0; i < vertex_cache.size(); ++i)
			{
				int new_id = vertex_cache[i]->id;
				int grow_id = grow_dir[new_id];
				if (!grow_flag[grow_id])
					continue;
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
			}

#if PRINT_DEBUG_INFO
			dprint("检测相似性能量");
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
					grow_flag[grow_dir[vc->id]] = false;
				}
			}
	    }
		else
			if_similarity_energy_low.resize(nvl, true);

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
					grow_flag[i] = false;
			}
		}

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
		//若已超过能量阈值或者遇到接近平面的区域，则退出
		if (!grow_flag[0] && !grow_flag[1])
			return false;

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
						return false;
					newv_flag[vvid] = true;
					new_vertex.push_back(&m4.verticelayers[vvid]);
					grow_dir[vvid] = growid;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}
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
							return false;
						newv_flag[vvid] = true;
						new_vertex.push_back(&m4.verticelayers[vvid]);
						grow_dir[vvid] = growid;
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);
			}
			begin_ = end_;
		}
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
							goto target2;
						ht = ht->next;
					} while (ht != hb);
					newf_flag[vfid] = true;
					new_face.push_back(&m4.facelayers[vfid]);
				target2:;
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
		vl = m4.conj_vl(vl, 1);
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
			if (!visited[m4.conj_vl(&m4.verticelayers[hl_mark->to], 3)->id])
				break;
			cut.push_back(m4.conj_vl(&m4.verticelayers[hl_mark->to], 3));
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
			if (!visited[m4.conj_vl(&m4.verticelayers[hl_mark->to], 3)->id])
				break;
			cut.push_back(m4.conj_vl(&m4.verticelayers[hl_mark->to], 3));
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

	void LoopGen::OptimizeLoop()
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

			LocalParametrization lp(m4, ip.vl);
			ConstructInitialRegion(&InfoOnMesh[ip.vl->id], lp);
			bool grow_flag[2] = { true, true };
			do
			{
				BoolVector visited_v = lp.region_v_flag;
				for (auto& ver : lp.new_vertex)
					visited_v[ver->id] = true;
				if (!ConstructRegionCut(ip.vl, visited_v, lp.cut))
					break;
				lp.run(cf->getNormal());
			} while (SpreadSubRegion(lp, grow_flag));

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

			if (lp.region_vertex.size() > 1)
			{
				seed_vertex.push_back(ip.vl);
				dprint("seed vertex:", seed_index, ip.vl->v.idx(), ip.vl->id, ip.data);
			}
			else
			{
				//dprint("undrafted vertex:", ip.vl->v.idx(), ip.vl->id, ip.energy);
			}
			if (lp.region_vertex.size() > 1)
			{
				for (auto rv : lp.region_vertex)
				{
					region_index[rv->id].push_back(seed_index);
					region_index[m4.conj_vl(rv, 2)->id].push_back(seed_index);
				}
			}
			for (auto rf : lp.region_face)
			{
				constraint_flag[rf->f.idx()] = true;
				constraint_dir.col(rf->f.idx()) = lp.x_axis.col(rf->f.idx());
			}
			auto& region_face = lp.region_face;
			auto& apl = lp.all_pl;
			for (int i = 0; i < apl.size(); ++i)
			{
				if (apl[i].empty())
					continue;
				//all_plane_loop[i] = std::move(apl[i]);
			}
#if 1
			if (lp.region_vertex.size() > 1)
			{
				cset.cylinders.push_back(cylinder());
				cylinder &cy = cset.cylinders.back();
				cy.id = cset.cylinders.size() - 1;
				cy.vertices = std::move(lp.region_vertex);
				cy.vertice_flag = std::move(lp.region_v_flag);
				/*cy.faces = std::move(lp.region_face);
				cy.cut = std::move(lp.cut);
				cy.vertice_flag = std::move(lp.region_v_flag);
				cy.face_flag = std::move(lp.region_f_flag);
				cy.cutv_flag = std::move(lp.cutv_flag);
				cy.cutf_flag = std::move(lp.cutf_flag);
				cy.vidmap = std::move(lp.vidmap);
				cy.uv[0] = std::move(lp.uv[0]);
				cy.uv[1] = std::move(lp.uv[1]);*/
				//cy.set_bound(lp.cut, lp.region_f_flag);
			}
#endif
			//break;
		}
		//return;
		ProcessOverlap(region_index);
		WriteRegion(cset.cylinders, cf->crossfield, model_name);

		cf->setOuterConstraint(constraint_flag, constraint_dir);
		cf->setField();
		cf->initFieldInfo();
		cf->write_field();
		tr.out("repair field:");
		dprint("complete");
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
		//int vvvid = InfoOnMesh[16389 * 2].energy < InfoOnMesh[16389 * 2 + 1].energy ? 16389 * 2 : 16389 * 2 + 1;
		for (int i = 0; i < mesh->n_vertices(); ++i)
		{
			if (m4.sing_flag[i])
				continue;
			for (int j = 0; j < 2; ++j)
			{
				//VertexLayer* vl = &m4.verticelayers[m4.verticemap[i] + j];
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

	bool LoopGen::CheckTopology(std::vector<VertexLayer*>& vertex_set, BoolVector& vs_flag, std::vector<int>& grow_dir)
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

	void LoopGen::ProcessOverlap(std::vector<std::vector<int>> &region_index)
	{
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
						oppo_vl = m4.conj_vl(vl, 2);
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
		int size = vh_set.size();
		cset.cylinders.swap(std::vector<cylinder>());
		cset.cylinders.resize(size);
		double halfSquare2 = sqrt(2)*0.5;
		cset.bound_edge_flag.resize(mesh->n_edges(), false);
		cset.vertex_bound_index.resize(mesh->n_vertices(), std::pair<int, int>(-1, -1));
		cset.has_vertex.resize(mesh->n_vertices(), false);
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
				seed_vl = m4.conj_vl(seed_vl, 1);

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
						if (vh_flag[vl->v.idx()] && !vl_flag[vl->id])
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
			ConstructRegionCut(seed_vl, vl_flag, cset.cylinders[i].cut);

			//初始化bound
			cset.cylinders[i].set_bound();

			//计算参数化
			cset.cylinders[i].parametrize(m4, cf->getNormal());

			//初始化index记录
			for (auto vl : vls)
			{
				cset.has_vertex[vl->v.idx()] = true;
			}
			for (int j = 0; j < 2; ++j)
			{
				for (auto hl : cset.cylinders[i].bounds[j])
				{
					cset.bound_edge_flag[hl->id / 8] = true;
					cset.vertex_bound_index[hl->to / 4].first = i;
					cset.vertex_bound_index[hl->to / 4].second = j;
				}
			}
		}

		//初始化各区域的交集情况
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
#if 1
		//从柱体边界出发搜索路径
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

		vl_pair_pq path_pq;
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				auto &bound = cset.cylinders[i].bounds[j];
				int shift = j == 0 ? 3 : 1;
				for (auto hl : bound)
				{
					VertexLayer* vl = m4.conj_vl(&m4.verticelayers[hl->to], shift);
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
					if (fromb != tob && !cset.tangential_intersection[fromb.first][tob.first])
						path_pq.emplace(vl, l);
				}
			}
		}
#endif
		//return;
		OptimizeDisk(path_pq);
	}

#if 1
	void LoopGen::OptimizeDisk(vl_pair_pq &path_pq)
	{
		int ccs = cset.cylinders.size();
		std::vector<BoolVector> flag(ccs);
		for (int i = 0; i < ccs; ++i)
			flag[i].resize(ccs, false);

		std::pair<int, int> from_bound;
		std::pair<int, int> to_bound;

		BoolVector constraint_flag(mesh->n_faces(), false);
		Eigen::Matrix3Xd constraint_dir(3, mesh->n_faces());
		all_path.resize(mesh->n_vertices());
		while (true)
		{
			vl_pair ip;
			ip = path_pq.top(); path_pq.pop();
			from_bound = cset.vertex_bound_index[ip.vl->v.idx()];
			to_bound = cset.vertex_bound_index[cset.all_path[ip.vl->id].back().hl->to / 4];
			while (flag[from_bound.first][to_bound.first])
			{
				if (path_pq.empty())
					goto goto0;
				ip = path_pq.top(); path_pq.pop();
				from_bound = cset.vertex_bound_index[ip.vl->v.idx()];
				to_bound = cset.vertex_bound_index[cset.all_path[ip.vl->id].back().hl->to / 4];
			}
			if (ip.vl->id == 119321)
				break;
			flag[from_bound.first][to_bound.first] = true;
			dprint(ip.vl->id, path_pq.size(), from_bound.first, to_bound.first);
			static double min_l = ip.data;
			if (ip.data > min_l * 5.0)
				break;

			bool grow_flag[2] = { true,true };
			temp_name tn;
			tn.from_bound = from_bound;
			tn.to_bound = to_bound;

			ConstructInitialRegion(ip.vl, tn);
			set_one_field(tn);
			while (true)
			{
				parame(tn);
				if (!SpreadSubRegion(tn, grow_flag))
					break;
				//break;
			}

			for (auto f : tn.region_face)
			{
				constraint_flag[f.idx()] = true;
				constraint_dir.col(f.idx()) = tn.x_axis.col(f.idx());
			}
			new_vert_flag = tn.new_v_flag;
			old_vert_flag = tn.region_v_flag;
			uv_para[0] = tn.uv[0];
			uv_para[1] = tn.uv[1];
			vertexidmap = tn.vidmap;
			new_face_flag = tn.new_f_flag;
			old_face_flag = tn.region_f_flag;
			xaxis = tn.x_axis;
			/*bool grow_flag[2] = { true, true };
			do
			{
				BoolVector visited_v = lp.region_v_flag;
				for (auto& ver : lp.new_vertex)
					visited_v[ver->id] = true;
				if (!ConstructRegionCut(ip.vl, visited_v, lp.cut))
					break;
				lp.run(cf->getNormal());
			} while (SpreadSubRegion(grow_flag));*/
			static int se = 0;
			++se;
			if (se > 4)
				break;
		}
	goto0:;
		return;
		for (auto &cy : cset.cylinders)
		{
			for (auto fh : cy.faces)
			{
				constraint_flag[fh->f.idx()] = true;
				constraint_dir.col(fh->f.idx()) = cf->crossfield.col(fh->id);
			}
		}
		cf->setOuterConstraint(constraint_flag, constraint_dir);
		cf->setField();
		cf->initFieldInfo();
		//cf->write_field();
		tr.out("repair field:");
		dprint("complete");
	}

	void LoopGen::ConstructInitialRegion(VertexLayer* vl, temp_name &tn)
	{
		auto &pl = cset.all_path[vl->id];
		auto &grow_dir = tn.grow_dir;
		grow_dir.resize(mesh->n_vertices(), -1);
		grow_dir[vl->v.idx()] = 0;
		/*dprint("vl:", vl->id / 4);
		for (auto hh : pl)
			dprint(hh.hl->id / 8);*/

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
		advancing_front[0].push_back(std::move(hierarchy_vertex[0]));
		advancing_front[1].push_back(std::move(hierarchy_vertex[1]));

		for (int itertimes = 0; itertimes < 3; ++itertimes)
		{
			const auto& current_af = advancing_front[0].back();
			std::vector<int> hierarchy;
			for (auto caf : current_af)
			{
				for (auto vv : mesh->vv_range(mesh->vertex_handle(caf)))
				{
					if (visited_vl[vv.idx()])
						continue;
					if (cset.has_vertex[vv.idx()] && cset.vertex_bound_index[vv.idx()].first < 0)
						continue;
					hierarchy.push_back(vv.idx());
					grow_dir[vv.idx()] = 0;
					visited_vl[vv.idx()] = true;
				}
			}
			advancing_front[0].push_back(std::move(hierarchy));
		}

		for (int itertimes = 0; itertimes < 3; ++itertimes)
		{
			const auto& current_af = advancing_front[1].back();
			std::vector<int> hierarchy;
			for (auto caf : current_af)
			{
				for (auto vv : mesh->vv_range(mesh->vertex_handle(caf)))
				{
					if (visited_vl[vv.idx()])
						continue;
					if (cset.has_vertex[vv.idx()] && cset.vertex_bound_index[vv.idx()].first < 0)
						continue;
					hierarchy.push_back(vv.idx());
					grow_dir[vv.idx()] = 1;
					visited_vl[vv.idx()] = true;
				}
			}
			advancing_front[1].push_back(std::move(hierarchy));
		}
		
		auto &new_vertex = tn.new_vertex;
		auto &newv_flag = tn.new_v_flag;
		auto &region_vertex = tn.region_vertex;
		auto &regionv_flag = tn.region_v_flag;
		//new_vertex.push_back(vl->v);
		//newv_flag.resize(nv, false); newv_flag[vl->v.idx()] = true;
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

		auto &new_face = tn.new_face;
		auto &newf_flag = tn.new_f_flag;
		newf_flag.resize(nf, false);
		for (auto new_v : new_vertex)
		{
			for (auto vh : mesh->voh_range(new_v))
			{
				if ((newv_flag[vh.to().idx()] || regionv_flag[vh.to().idx()]) &&
					(newv_flag[vh.next().to().idx()] || regionv_flag[vh.next().to().idx()]))
				{
					if (!newf_flag[vh.face().idx()])
					{
						new_face.push_back(vh.face());
						newf_flag[vh.face().idx()] = true;
					}
				}
			}
		}

		//加入边界上x_axis, y_axis的约束
		auto &x_axis = tn.x_axis;
		auto &y_axis = tn.y_axis;
		x_axis.resize(3, mesh->n_faces()); x_axis.setZero();
		y_axis.resize(3, mesh->n_faces()); y_axis.setZero();
		auto &constraint_f_flag = tn.constraint_f_flag;
		constraint_f_flag.resize(mesh->n_faces(), false);
		auto &crossfield = cf->crossfield;
		//可参考《实验记录》2.25示例图
		auto &bound0 = cset.cylinders[tn.from_bound.first].bounds[tn.from_bound.second];
		for (auto hl : bound0)
		{
			int flid = hl->oppo->left;
			if (tn.from_bound.second == 1)
				flid = (flid / 4) * 4 + ((flid + 2) % 4);
			y_axis.col(flid / 4) = crossfield.col(flid);
			x_axis.col(flid / 4) = crossfield.col((flid / 4) * 4 + ((flid + 3) % 4));
			constraint_f_flag[flid / 4] = true;
		}
		auto &bound1 = cset.cylinders[tn.to_bound.first].bounds[tn.to_bound.second];
		for (auto hl : bound1)
		{
			int flid = hl->oppo->left;
			if (tn.to_bound.second == 0)
				flid = (flid / 4) * 4 + ((flid + 2) % 4);
			y_axis.col(flid / 4) = crossfield.col(flid);
			x_axis.col(flid / 4) = crossfield.col((flid / 4) * 4 + ((flid + 3) % 4));
			constraint_f_flag[flid / 4] = true;
		}

		tn.vidmap.resize(nv, -1);
		tn.vidmap[vl->v.idx()] = 0;
		tn.uv[0].resize(1); tn.uv[0](0) = 0;
		tn.uv[1].resize(1); tn.uv[1](0) = 0;
		tn.region_f_flag.resize(mesh->n_faces(), false);
	}

	void LoopGen::set_one_field(temp_name &tn)
	{
		auto &new_vertex = tn.new_vertex;
		auto &newv_flag = tn.new_v_flag;
		auto &new_face = tn.new_face;
		auto &newf_flag = tn.new_f_flag;
		auto &x_axis = tn.x_axis;
		auto &y_axis = tn.y_axis;

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
				if (newf_flag[gid] && fid < gid || tn.constraint_f_flag[gid])
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

	void LoopGen::parame(temp_name &tn)
	{
		auto &new_vertex = tn.new_vertex;
		auto &newv_flag = tn.new_v_flag;
		auto &new_face = tn.new_face;
		auto &newf_flag = tn.new_f_flag;
		auto &x_axis = tn.x_axis;
		auto &y_axis = tn.y_axis;
		auto &region_vertex = tn.region_vertex;
		auto &regionv_flag = tn.region_v_flag;
		auto &region_face = tn.region_face;
		auto &regionf_flag = tn.region_f_flag;
		auto &vidmap = tn.vidmap;
		Eigen::VectorXd* uv = tn.uv;

		std::vector<std::map<VertexHandle, Eigen::Vector3d>> info(mesh->n_faces());
		auto &normal = cf->getNormal();
		HalfedgeHandle hh;
		FaceHandle fh;
		int fid;
		auto calc_vector = [&]()
		{
			double inv_area = 1.0 / (2 * mesh->calc_face_area(fh));
			auto vi = mesh->to_vertex_handle(hh);
			auto ev = mesh->calc_edge_vector(mesh->prev_halfedge_handle(hh));
			info[fid].insert(std::make_pair(vi, normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area));
			vi = mesh->to_vertex_handle(mesh->next_halfedge_handle(hh));
			ev = mesh->calc_edge_vector(hh);
			info[fid].insert(std::make_pair(vi, normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area));
			vi = mesh->from_vertex_handle(hh);
			ev = mesh->calc_edge_vector(mesh->next_halfedge_handle(hh));
			info[fid].insert(std::make_pair(vi, normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area));
		};
		for (auto f : new_face)
		{
			hh = mesh->fh_begin(f).handle();
			fh = mesh->face_handle(hh);
			fid = fh.idx();
			calc_vector();
		}

		int region_vertex_size = region_vertex.size();
		//计算idmap
		{
			int count = region_vertex_size;
			for (auto v : new_vertex)
			{
				vidmap[v.idx()] = count++;
			}
		}

		int new_vertex_size = new_vertex.size();
		int new_face_size = new_face.size();
		std::vector<Eigen::Triplet<double>> triple;
		std::vector<double> w(mesh->n_vertices(), 0);
		//std::vector<double> size_ratio(nf, 1.0);
		uv[0].conservativeResize(region_vertex_size + new_vertex_size); uv[0].tail(new_vertex_size).setZero();
		uv[1].conservativeResize(region_vertex_size + new_vertex_size); uv[1].tail(new_vertex_size).setZero();

		for (auto v : new_vertex)
		{
			int vid = v.idx();
			int vm = vidmap[vid];
			for (auto vf : mesh->vf_range(v))
			{
				int vf_id = vf.idx();
				if (!newf_flag[vf_id])
					continue;
				Eigen::Vector3d &R = info[vf_id][v];
				for (const auto &f_info : info[vf_id])
				{
					double dot_ = R.dot(f_info.second);
					int fvid = f_info.first.idx();
					if (!newv_flag[fvid])
					{
						uv[0](vm) -= dot_ * uv[0](vidmap[fvid]);
						uv[1](vm) -= dot_ * uv[1](vidmap[fvid]);
					}
					else
					{
						w[fvid] += dot_;
					}
					uv[0](vm) += x_axis.col(vf_id).dot(R);
					uv[1](vm) += y_axis.col(vf_id).dot(R);
				}
			}

			triple.emplace_back(vm - region_vertex_size, vm - region_vertex_size, w[vid]);
			w[vid] = 0;
			for (auto vv : mesh->vv_range(v))
			{
				int vvid = vv.idx();
				if (!newv_flag[vvid])
					continue;
				triple.emplace_back(vm - region_vertex_size, vidmap[vvid] - region_vertex_size, w[vvid]);
				w[vvid] = 0;
			}
		}

		Eigen::SparseMatrix<double> A(new_vertex_size, new_vertex_size);
		A.setFromTriplets(triple.begin(), triple.end());
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		uv[0].tail(new_vertex_size) = solver.solve(uv[0].tail(new_vertex_size));
		uv[1].tail(new_vertex_size) = solver.solve(uv[1].tail(new_vertex_size));
	}

	bool LoopGen::RefinePathByParametrization(VertexHandle v, temp_name &tn, BoolVector &visited_v, BoolVector &visited_f)
	{
		double v_para = tn.GetV(v.idx());
		VertexLayer* vl = &m4.verticelayers[v.idx() * 4];

		auto hl_begin = vl->hl;
		auto hl_transfer = hl_begin;
		while (!visited_f[hl_transfer->left / 4])
		{
			hl_transfer = hl_transfer->prev->oppo;
		}
		hl_begin = hl_transfer;
		double s0 = tn.GetV(hl_transfer->to / 4) - v_para;
		PointOnHalfedgeLayer pohl[2];
		int id[2] = { 0,0 };

		double s1 = tn.GetV(hl_transfer->next->to / 4) - v_para;
		do
		{
			if (visited_v[hl_transfer->next->to / 4])
				s1 = tn.GetV(hl_transfer->next->to / 4) - v_para;
			if (visited_f[hl_transfer->left / 4])
			{
				if (s0*s1 < 0)
				{
					if (s0 > 0)
					{
						pohl[0].hl = hl_transfer->next->oppo;
						pohl[0].c = s0 / (s0 - s1);
						++id[0];
					}
					else
					{
						pohl[1].hl = hl_transfer->next;// ->oppo;
						pohl[1].c = s1 / (s1 - s0);
						++id[1];
					}
				}
			}
			s0 = s1;
			hl_transfer = hl_transfer->prev->oppo;
		} while (hl_transfer->id / 4 != hl_begin->id / 4);

		if (id[0] > 1 || id[1] > 1)
			return false;
		if (id[0] != 1 && cset.vertex_bound_index[v.idx()] != tn.from_bound)
			return false;
		if (id[1] != 1 && cset.vertex_bound_index[v.idx()] != tn.to_bound)
			return false;

		PlaneLoop planar_path;
		double distance[2] = { 0,0 };
		if (cset.vertex_bound_index[v.idx()] != tn.from_bound)
		{
			planar_path.push_back(pohl[0]);
			distance[0] = tn.GetV(pohl[0].hl->from / 4) - v_para;
			distance[1] = tn.GetV(pohl[0].hl->to / 4) - v_para;

			hl_transfer = pohl[0].hl;
			while (!cset.bound_edge_flag[hl_transfer->id / 8])
			{
				int vid = hl_transfer->next->to / 4;
				if (!visited_v[vid])
					return false;
				double s = tn.GetV(vid) - v_para;
				if (s > 0)
				{
					hl_transfer = hl_transfer->prev->oppo;
					distance[1] = s;
				}
				else
				{
					hl_transfer = hl_transfer->next->oppo;
					distance[0] = s;
				}
				planar_path.emplace_back(hl_transfer, distance[1] / (distance[1] - distance[0]));
			}
		}
		std::reverse(planar_path.begin(), planar_path.end());
		planar_path.emplace_back(nullptr, v.idx() * 1.0);
		if (cset.vertex_bound_index[v.idx()] != tn.to_bound)
		{
			planar_path.push_back(pohl[1]);
			distance[0] = tn.GetV(pohl[1].hl->from / 4) - v_para;
			distance[1] = tn.GetV(pohl[1].hl->to / 4) - v_para;

			hl_transfer = pohl[1].hl;
			while (!cset.bound_edge_flag[hl_transfer->id / 8])
			{
				int vid = hl_transfer->oppo->next->to / 4;
				if (!visited_v[vid])
					return false;
				double s = tn.GetV(vid) - v_para;
				if (s > 0)
				{
					hl_transfer = hl_transfer->oppo->next;
					distance[1] = s;
				}
				else
				{
					hl_transfer = hl_transfer->oppo->prev;
					distance[0] = s;
				}
				planar_path.emplace_back(hl_transfer, distance[1] / (distance[1] - distance[0]));
			}
		}
		all_path[v.idx()] = std::move(planar_path);
		//lp.all_pl[vl->id] = std::move(planar_loop);
		return true;
	}

	double LoopGen::AssembleSimilarityAngle(VertexHandle v, Eigen::VectorXd& sa, temp_name &tn, int path_fragment_num)
	{
		Eigen::Matrix3Xd fragment(3, path_fragment_num); fragment.setZero();
		PointOnHalfedgeLayer pohl;
		double fromu = 0, tou = 0;
		OpenMesh::Vec3d frompos, topos, ev;
		double step = path_fragment_num / (tn.umx[1] - tn.umx[0]);
		double len = 0;
		auto setData = [&]()
		{
			if (tou <= tn.umx[0] || fromu >= tn.umx[1]) return;
			if (fromu < tn.umx[0]) fromu = tn.umx[0];
			if (tou > tn.umx[1]) tou = tn.umx[1];
			int u0 = std::floor((fromu - tn.umx[0]) * step);
			int u1 = std::floor((tou - tn.umx[0]) * step);
			if (u1 >= path_fragment_num) u1 = path_fragment_num - 1;
			if (u0 == u1) return;
			ev = topos - frompos;
			len += ev.norm();
			ev.normalized();
			for (int i = u0 + 1; i <= u1; ++i) fragment.col(i) << ev[0], ev[1], ev[2];
		};
		frompos = all_path[v.idx()].front().point(m4);
		pohl = all_path[v.idx()].front();
		if (pohl.hl) fromu = pohl.c*tn.GetU(pohl.hl->from / 4) + (1 - pohl.c)*tn.GetU(pohl.hl->to / 4);
		else fromu = tn.GetU(pohl.c);
		for (int i = 1; i < all_path[v.idx()].size(); ++i)
		{
			pohl = all_path[v.idx()][i];
			topos = pohl.point(m4);
			if (pohl.hl) tou = pohl.c*tn.GetU(pohl.hl->from / 4) + (1 - pohl.c)*tn.GetU(pohl.hl->to / 4);
			else tou = tn.GetU(pohl.c);
			setData();
			frompos = topos;
			fromu = tou;
		}
		sa.resize(path_fragment_num * (path_fragment_num + 1) / 2); sa.setZero();
		int count = 0;
		for (int i = 0; i < path_fragment_num; ++i)
			for (int j = i + 1; j < path_fragment_num; ++j)
				sa(count++) = conservativeArcCos(fragment.col(i).dot(fragment.col(j)));
		return len;
	}

	bool LoopGen::SpreadSubRegion(temp_name &tn, bool grow_flag[2])
	{
		int nv = mesh->n_vertices();
		int nf = mesh->n_faces();

		auto& new_vertex = tn.new_vertex;
		auto& new_face = tn.new_face;
		auto &region_vertex = tn.region_vertex;
		auto &region_face = tn.region_face;
		auto& newv_flag = tn.new_v_flag;
		auto& newf_flag = tn.new_f_flag;
		auto &regionv_flag = tn.region_v_flag;
		auto &regionf_flag = tn.region_f_flag;
		BoolVector visited_v = tn.region_v_flag;
		BoolVector visited_f = tn.region_f_flag;
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
			RefinePathByParametrization(region_vertex.front(), tn, visited_v, visited_f);
		std::vector<VertexHandle> vertex_cache;
		BoolVector vertex_cache_flag = tn.region_v_flag;
		//return true;
		for (auto new_v : new_vertex)
		{
			//int se = tn.all_pl[new_v->id].size();
			if (RefinePathByParametrization(new_v, tn, visited_v, visited_f))
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
		if (vertex_cache.empty())
			return false;

#if PRINT_DEBUG_INFO
		dprint("拓扑检查");
#endif

		auto& grow_dir = tn.grow_dir;

		//检测相似性能量
		auto& normal_similarity_angle = tn.normal_similarity_angle;
		int path_fragment_num = all_path[region_vertex.front().idx()].size();
		//auto& all_pl = lp.all_pl;
		//int loop_fragment_num = all_pl[lp.region_vertex.front()->id].size();
		if (vertex_cache_flag[region_vertex.front().idx()] && !tn.has_nsa)
		{
			tn.has_nsa = true;
			tn.umx[0] = tn.GetU(region_vertex.front().idx());
			auto pohl = all_path[region_vertex.front().idx()].back();
			tn.umx[1] = pohl.c*tn.GetU(pohl.hl->from / 4) + (1 - pohl.c)*tn.GetU(pohl.hl->to / 4);
			double len = AssembleSimilarityAngle(region_vertex.front(), normal_similarity_angle, tn, path_fragment_num);
			tn.length[0] = std::min(tn.length[0], len);
			tn.length[1] = std::max(tn.length[1], len);
			tn.length[2] = std::max(tn.length[2], len);
		}

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
				Eigen::VectorXd similarity_angle;
				double len = AssembleSimilarityAngle(vertex_cache[i], similarity_angle, tn, path_fragment_num);
				double sum = 0;
				double seaa = 0;
				for (int j = 0; j < similarity_angle.size(); ++j) {
					sum += fabs(normal_similarity_angle(j) - similarity_angle(j));
					seaa = std::max(seaa, fabs(normal_similarity_angle(j) - similarity_angle(j)));
				}
				tn.length[0] = std::min(tn.length[0], len);
				tn.length[grow_id + 1] = std::max(tn.length[grow_id + 1], len);
				if (sum < disk_e * similarity_angle.size() && tn.length[0] * 1.5 > tn.length[grow_id + 1])
				{
					if_similarity_energy_low[new_id] = true;
				}
				dprint(i, sum / similarity_angle.size(), seaa, disk_e);
			}

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
			//if_similarity_energy_low.resize(nv, true);

		//for (int i = 0; i < 2; ++i)
		//{
		//	if (grow_flag[i])
		//	{
		//		double grad = LoopLenGrad(vertex_cache, lp, vertex_cache_flag, i);
		//		//dprint("grad:", i, grad);
		//		if (grad > PI)
		//			grow_flag[i] = false;
		//	}
		//}

		//更新新区域
		int count = region_vertex.size();
		int begin_ = count;
		auto& u_para = tn.GetU();
		auto& v_para = tn.GetV();
		auto& vidmap = tn.vidmap;
		for (auto new_v : vertex_cache)
		{
			if (if_similarity_energy_low[new_v.idx()])
			{
				region_vertex.push_back(new_v);
				int newid = new_v.idx();
				regionv_flag[newid] = true;
				u_para(count) = tn.GetU(newid);
				v_para(count) = tn.GetV(newid);
				vidmap[newid] = count++;
			}
		}
		u_para.conservativeResize(count);
		v_para.conservativeResize(count);

		for (auto new_f : new_face)
		{
			/*auto hl_begin = new_f->hl;
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
		target1:;*/
			for (auto fv : mesh->fv_range(new_f))
			{
				if (!regionv_flag[fv.idx()])
					goto target1;
			}
			regionf_flag[new_f.idx()] = true;
			region_face.push_back(new_f);
		target1:;
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
			/*auto hl_begin = rvi->hl;
			auto hl_transfer = hl_begin;
			do
			{
				int vvid = hl_transfer->to;
				if (!(regionv_flag[vvid] || newv_flag[vvid]))
				{
					if (m4.sing_flag[m4.verticelayers[vvid].v.idx()])
						return false;
					newv_flag[vvid] = true;
					new_vertex.push_back(&m4.verticelayers[vvid]);
					grow_dir[vvid] = growid;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);*/
			for (auto vv : mesh->vv_range(rvi))
			{
				int vvid = vv.idx();
				if (regionv_flag[vvid] || newv_flag[vvid])
					continue;
				if (cset.has_vertex[vvid] && cset.vertex_bound_index[vvid].first < 0)
					continue;
				new_vertex.push_back(vv);
				newv_flag[vvid] = true;
				grow_dir[vvid] = growid;
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
				/*auto hl_begin = rvj->hl;
				auto hl_transfer = hl_begin;
				do
				{
					int vvid = hl_transfer->to;
					if (!(regionv_flag[vvid] || newv_flag[vvid]))
					{
						if (m4.sing_flag[m4.verticelayers[vvid].v.idx()])
							return false;
						newv_flag[vvid] = true;
						new_vertex.push_back(&m4.verticelayers[vvid]);
						grow_dir[vvid] = growid;
					}
					hl_transfer = hl_transfer->prev->oppo;
				} while (hl_transfer != hl_begin);*/
				for (auto vv : mesh->vv_range(rvj))
				{
					int vvid = vv.idx();
					if (regionv_flag[vvid] || newv_flag[vvid])
						continue;
					if (cset.has_vertex[vvid] && cset.vertex_bound_index[vvid].first < 0)
						continue;
					new_vertex.push_back(vv);
					newv_flag[vvid] = true;
					grow_dir[vvid] = growid;
				}
			}
			begin_ = end_;
		}

		newf_flag.resize(nf, false);
		new_face.clear();
		for (auto new_v : new_vertex)
		{
			int newid = new_v.idx();
			for (auto voh : mesh->voh_range(new_v))
			{
				int vfid = voh.face().idx();
				if (regionf_flag[vfid] || newf_flag[vfid])
					continue;
				if (!regionv_flag[voh.to().idx()] && !newv_flag[voh.to().idx()])
					continue;
				if (!regionv_flag[voh.next().to().idx()] && !newv_flag[voh.next().to().idx()])
					continue;
				new_face.push_back(voh.face());
				newf_flag[voh.face().idx()] = true;
			}
			/*auto hl_begin = new_v->hl;
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
							goto target2;
						ht = ht->next;
					} while (ht != hb);
					newf_flag[vfid] = true;
					new_face.push_back(&m4.facelayers[vfid]);
				target2:;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);*/
		}
#if PRINT_DEBUG_INFO
		dprint("扩展advancing_front");
#endif

		//优化新区域的场
		for (auto rf : region_face)
			tn.constraint_f_flag[rf.idx()] = true;
		set_one_field(tn);
#if PRINT_DEBUG_INFO
		dprint("优化新区域的场");
#endif
#endif
		return true;
	}
#endif

}