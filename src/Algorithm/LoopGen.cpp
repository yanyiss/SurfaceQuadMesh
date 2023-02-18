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
#if 1
			fragment.col(0) = vec;
#else
			fragment.col(0) = loop.col(mark % cols) + u0 * vec;
#endif
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
#if 1
				fragment.col(i) = vec;
#else
				fragment.col(i) = loop.col((mark + 1) % cols) + (u0 - l) * vec;
#endif
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
				//dot = fabs(fragment0.col(i).dot(fragment0.col(j)))/(fabs((fragment1.col(i).dot(fragment1.col(j)))) + YYSS_FAIRLY_SMALL);
				//sum += fabs(sin(normal_similarity_angle(j) - similarity_angle(j))) + 1 - cos(normal_similarity_angle(j) - similarity_angle(j));
				double theta = conservativeArcCos(fragment0.col(i).dot(fragment0.col(j))) - conservativeArcCos(fragment1.col(i).dot(fragment1.col(j)));
				sum += fabs(theta);
				//sum += std::min(100.0, dot + 1.0 / (dot + YYSS_FAIRLY_SMALL) - 2);
				//dprint(i, j, fabs(fragment0.col(i).dot(fragment0.col(j))), fabs((fragment1.col(i).dot(fragment1.col(j)))), dot, dot + 1.0 / (dot + YYSS_FAIRLY_SMALL) - 2);
				//dot = fragment0.col(i).dot(fragment0.col(j)) - fragment1.col(i).dot(fragment1.col(j));
				//sum += dot * dot;
			}
		}
		return 2.0 * sum / (n * (n - 1));
		//return 0;
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
		//double length = 0;
		for (int i = 0; i < n; ++i)
		{
			sum += fabs(xyz[0](i) * plane[0] + xyz[1](i) * plane[1] + xyz[2](i) * plane[2] + plane[3]);
			//length += (position.col(i) - position.col((i + 1) % n)).norm();
		}
		//dprint("eee:", sum / n, sum, n, length);
		//return sum / (n * length);
		return sum / n;
	}

	double LoopGen::RefineLoopByPlanarity(/*M4 &m4, */std::vector<VertexLayer*>& loop, PlaneLoop& planar_loop)
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
			//dprint(m4.verticelayers[hl_transfer->to].v.idx());
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

		/*for (; hitr != mesh->voh_end(loop[0]); ++hitr)
		{
			auto s1 = dis(mesh->point(mesh->to_vertex_handle(hitr.handle())));
			if (s0 * s1 < 0)
			{
				if (s0 > 0)
				{
					poh[0].h = mesh->next_halfedge_handle(hitr.handle());
					poh[0].c = s0 / (s0 - s1);
					distance[0] = s0; distance[1] = s1;
					++id[0];
				}
				else
				{
					poh[1].h = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(hitr.handle()));
					poh[1].c = s1 / (s1 - s0);
					++id[1];
				}
			}
			s0 = s1;
		}*/
		if (id[0] != 1 || id[1] != 1)
		{
			//dprint("error in repairing loop:", loop[0].idx());
			return YYSS_INFINITE;
		}

		//检查出发点到pohl[0]的方向与搜索loop的方向是否相同，若不相同，则调换pohl[0]和pohl[1]
		auto &v0 = cf->getCrossField().col(pohl[0].hl->left);
		//auto v1 = pohl[0].c*mesh->point(m4.verticelayers[pohl[0].hl->from].v)
		//	+ (1 - pohl[0].c)*mesh->point(m4.verticelayers[pohl[1].hl->to].v);
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
		//{
		//	auto& matching = cf->getMatching();
		//	//int shift = 0;
		//	auto h_end = mesh->prev_halfedge_handle(poh[0].h);
		//	auto h_begin = mesh->voh_begin(mesh->from_vertex_handle(h_end)).handle();
		//	while (h_begin != h_end)
		//	{
		//		h_begin = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h_begin));
		//		shift += 4 - matching[h_begin.idx()];
		//	}
		//	auto& v0 = cf->getCrossField().col(4 * mesh->face_handle(h_end).idx() + (shift % 4));
		//	auto v1 = poh[0].c * mesh->point(mesh->from_vertex_handle(poh[0].h)) +
		//		(1 - poh[0].c) * mesh->point(mesh->to_vertex_handle(poh[0].h)) - mesh->point(mesh->from_vertex_handle(h_end));
		//	if (v0(0) * v1[0] + v0(1) * v1[1] + v0(2) * v1[2] < 0)
		//	{
		//		std::swap(poh[0], poh[1]);
		//		for (int i = 0; i < 4; ++i)
		//			plane[i] *= -1.0;
		//		for (int i = 0; i < 2; ++i)
		//		{
		//			poh[i].h = mesh->opposite_halfedge_handle(poh[i].h);
		//			poh[i].c = 1 - poh[i].c;
		//		}
		//		distance[0] = dis(mesh->point(mesh->to_vertex_handle(poh[0].h)));
		//		distance[1] = dis(mesh->point(mesh->from_vertex_handle(poh[0].h)));
		//	}
		//}

		//planar_loop.push_back(poh[0]);
		//auto h = poh[0].h;
		//while (h.idx() != poh[1].h.idx())
		//{
		//	auto s = dis(mesh->point(mesh->from_vertex_handle(mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(h)))));
		//	if (s > 0)
		//	{
		//		h = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h));
		//		distance[0] = s;
		//	}
		//	else
		//	{
		//		h = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(h));
		//		distance[1] = s;
		//	}
		//	planar_loop.emplace_back(h, distance[0] / (distance[0] - distance[1]));
		//	//dprint(h.idx(),h.idx()/2);
		//	/*if (h.idx() == 17413)
		//		break;*/
		//}

		
		//planar_loop.push_back(pohl[0]);
		planar_loop.emplace_back(pohl[0].hl, pohl[0].c);
		hl_transfer = pohl[0].hl;
		//while (hl_transfer != pohl[1].hl)
		while (hl_transfer->h != pohl[1].hl->h)
		{
			//dprint(hl_transfer->h.idx() / 2);
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
			//planar_loop.emplace_back(hl_transfer, distance[0] / (distance[0] - distance[1]));
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
		//tr.begin("1");
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
					vid = m4.verticemap[vid];
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
		//tr.end("1");
		//tr.begin("2");
		//int siseaae = 0;
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

			//++siseaae;
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
							//dprint(pq.size());
							prev[toid] = hl_transfer;
							//visited[vid] = true;
							vid = m4.verticemap[vid];
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
		//dprint("siseaae:", siseaae);
		//tr.end("2");
		//tr.begin("3");
		loop.clear();
		loop.push_back(vl);
		VertexLayer* prev_layer = &m4.verticelayers[prev[vl->id]->from];
		while (prev_layer != vl)
		{
			loop.push_back(prev_layer);
			prev_layer = &m4.verticelayers[prev[prev_layer->id]->from];
		}
		loop.push_back(vl);
		//tr.end("3");
		//tr.begin("4");
		//dprint(loop.size());
		std::reverse(loop.begin(), loop.end());
		//tr.end("4");
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
		m4.update(); 
		m4.set_weight();
		//InfoOnMesh.resize(mesh->n_vertices() * 2);
		InfoOnMesh.resize(m4.verticelayers.size());
		int count = 0;
		for (int i = 0; i < InfoOnMesh.size(); ++i)
		{
			InfoOnMesh[i].id = i;
			InfoOnMesh[i].energy = YYSS_INFINITE;
			/*InfoOnMesh[i].plid = count;
			if (m4.verticelayers[i].layer % 2 == 0)
			{
				InfoOnMesh[i].dir = Forward;
				if (m4.sing_flag[m4.verticelayers[i].v.idx()])
					++count;
			}
			else
			{
				InfoOnMesh[i].dir = Reverse;
				++count;
			}*/
			if (m4.sing_flag[m4.verticelayers[i].v.idx()])
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
			}
		}
		int nv = mesh->n_vertices();

		pls.resize((m4.verticelayers.size() - cf->getSingularity().size()) / 2 + cf->getSingularity().size());

		tr.tog();
		//if (COMPUTE_NEW_PLANELOOP || !ReadPlaneLoop(m4, InfoOnMesh, model_name, mesh))
		if (COMPUTE_NEW_PLANELOOP || !ReadPlaneLoop(m4, pls, model_name, mesh))
		{
#pragma omp parallel for
			for (int i = 0; i < nv; ++i)
			{
				if (m4.sing_flag[i])
					continue;
				VertexLayer* vl = &m4.verticelayers[m4.verticemap[i]];
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
					//auto pos = mesh->point(mesh->to_vertex_handle(itr->h)) * (1 - itr->c) + mesh->point(mesh->from_vertex_handle(itr->h)) * itr->c;
					auto pos = itr->point(m4);
					loop.col(++c) << pos[0], pos[1], pos[2];
				}
			}
			else
			{
				for (auto itr = pl.rbegin(); itr != pl.rend(); ++itr)
				{
					//auto pos = mesh->point(mesh->to_vertex_handle(itr->h)) * (1 - itr->c) + mesh->point(mesh->from_vertex_handle(itr->h)) * itr->c;
					auto pos = itr->point(m4);
					loop.col(++c) << pos[0], pos[1], pos[2];
				}
			}
		};

		/*for (int i = 0; i < m4.xyvertices.size(); ++i)
		{
			VertexLayer* vl = m4.xyvertices[i];
			if (m4.sing_flag[vl->v.idx()])
			{
				HalfedgeLayer* hl = vl->hl;
				HalfedgeLayer* hl_transfer = hl;
				do
				{
					int toid = hl_transfer->to;
					if (vl->id > toid)
					{

					}
					hl_transfer = hl_transfer->prev->next;
				} while (hl_transfer != hl);
			}
		}*/
		//return;
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
					//similarity_energy[hl.oppo->id + ((hl.oppo->id + 2) % 4)] = similarity_energy[hlid];
					similarity_energy[(hl.oppo->id / 4) * 4 + ((hl.oppo->id + 2) % 4)] = similarity_energy[hlid];
					//dprint(similarity_energy[hl.id]);
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
		//advancing_front[0].push_back(std::vector<int>(1, iov->id));
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
		//auto& grow_dir = lp.GetGrowDir();
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
						//ht = ht->prev->oppo;
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

		/*static int rrr = 0;
		++rrr;
		if (rrr == 6)
		{
			int p = 0;
			for (auto new_v : lp.new_vertex)
				if (new_v->v.idx() == 20723)
					dprint("ueis", new_v->id);
		}*/
		//dprint("cut", lp.cut.size());
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
		////dprint("add vertex:", count - begin_);
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
		//dprint("find vertex:", new_vertex.size());
		//dprint("feayfusk", lp.new_v_flag[91260], lp.region_v_flag[91260]);
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

		/*if (rrr == 5)
		{
			int p = 0;
			return false;
		}*/
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
				if (!mesh->is_boundary(hl_transfer->h))
				{
					int gid = mesh->face_handle(hl_transfer->oppo->h).idx();
					int glid = hl_transfer->oppo->left;
					if (constraint_flag[glid] || (opt_flag[glid] && flid > glid))
					//if (!(!(opt_flag[glid] || constraint_flag[glid]) || (opt_flag[glid] && flid < glid)))
					{
						auto ev = (position.col(mesh->to_vertex_handle(hl_transfer->h).idx()) -
							position.col(mesh->from_vertex_handle(hl_transfer->h).idx())).normalized();
						COMPLEX e_f = COMPLEX(ev.dot(faceBase.col(fid * 2)), -ev.dot(faceBase.col(fid * 2 + 1)));
						//e_f /= std::sqrt(e_f.imag()*e_f.imag() + e_f.real()*e_f.real());
						COMPLEX e_g = COMPLEX(ev.dot(faceBase.col(gid * 2)), -ev.dot(faceBase.col(gid * 2 + 1)));
						//e_g /= std::sqrt(e_g.imag()*e_g.imag() + e_g.real()*e_g.real());
						if (opt_flag[flid])
						{
							triple.emplace_back(count, fidmap[flid], e_f);
						}
						else
						{
							COMPLEX dir = COMPLEX(x_axis.col(fid).dot(faceBase.col(2 * fid)), x_axis.col(fid).dot(faceBase.col(2 * fid + 1)));
							//COMPLEX dir = COMPLEX(crossfield.col(flid).dot(faceBase.col(fid * 2)), crossfield.col(flid).dot(faceBase.col(fid * 2 + 1)));
							b(count) -= e_f * dir;
						}
						if (opt_flag[glid])
						{
							triple.emplace_back(count, fidmap[glid], -e_g);
						}
						else
						{
							COMPLEX dir = COMPLEX(x_axis.col(gid).dot(faceBase.col(2 * gid)), x_axis.col(gid).dot(faceBase.col(2 * gid + 1)));
							//COMPLEX dir = COMPLEX(crossfield.col(glid).dot(faceBase.col(gid * 2)), crossfield.col(glid).dot(faceBase.col(gid * 2 + 1)));
							b(count) += e_g * dir;
						}
						++count;
					}
				}
				//hl_transfer = hl_transfer->prev->oppo;
				hl_transfer = hl_transfer->next;
			} while (hl_transfer != hl_begin);
		}
		//Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<COMPLEX>> solver;
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
			//dprint(hl_mark->to);
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
			//cut.push_back(&m4.verticelayers[hl_mark->to]);
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
		info_pair_pq pq;
		AssembleIOVLoopEnergy(pq);

		BoolVector ifset_flag(m4.verticelayers.size(), false);
		BoolVector constraint_flag(mesh->n_faces(), false);
		Eigen::Matrix3Xd constraint_dir(3, mesh->n_faces());
		all_plane_loop.resize(m4.verticelayers.size());

		while (true)
		{
			info_pair ip;
			ip = pq.top(); pq.pop();
			while (ifset_flag[ip.vl->id])
			{
				if (pq.empty())
				{
					ip.energy = YYSS_INFINITE;
					break;
				}
				ip = pq.top(); pq.pop();
			}
			if (ip.energy > energy_threshold)
				break;
			//ip.vl = &m4.verticelayers[151815];

			LocalParametrization lp(m4, ip.vl);
			ConstructInitialRegion(&InfoOnMesh[ip.vl->id], lp);

			bool grow_flag[2] = { true, true };
			int es = 0;
			do
			{
				++es;
				//dprint();
				//dprint(es, grow_flag[0], grow_flag[1]);
				BoolVector visited_v = lp.region_v_flag;
				for (auto& ver : lp.new_vertex)
					visited_v[ver->id] = true;

				if (!ConstructRegionCut(ip.vl, visited_v, lp.cut))
					break;
				//if (ip.vl->id == 151815)
				///	break;
				lp.run(cf->getNormal());
				//if (es > 1)
				//	break;
				//dprint("cut size:", lp.cut.size());
				for (int ii = 0; ii < lp.cut.size(); ++ii)
				{
					//dprint(ii, lp.cut[ii]->id);
				}
			} while (SpreadSubRegion(lp, grow_flag));

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

			if (lp.region_vertex.size() > 1)
			{
				seed_vertex.push_back(ip.vl);
				dprint("seed vertex:", seed_vertex.size() - 1, ip.vl->v.idx(), ip.vl->id, ip.energy);
			}
			else
			{
				dprint("undrafted vertex:", ip.vl->v.idx(), ip.vl->id, ip.energy);
			}
			for (auto rv : lp.region_vertex)
			{
				ifset_flag[rv->id] = true;
				ifset_flag[m4.conj_vl(rv, 2)->id] = true;
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
				all_plane_loop[i] = std::move(apl[i]);
			}
			/*if (lp.region_vertex.size() > 1)
			{
				cset.cylinders.push_back(cylinder());
				cylinder &cy = cset.cylinders.back();
				cy.id = cset.cylinders.size() - 1;
				cy.vertices = std::move(lp.region_vertex);
				cy.faces = std::move(lp.region_face);
				cy.cut = std::move(lp.cut);
				cy.vertice_flag = std::move(lp.region_v_flag);
				cy.face_flag = std::move(lp.region_f_flag);
				cy.cutv_flag = std::move(lp.cutv_flag);
				cy.cutf_flag = std::move(lp.cutf_flag);
				cy.vidmap = std::move(lp.vidmap);
				cy.uv[0] = std::move(lp.uv[0]);
				cy.uv[1] = std::move(lp.uv[1]);
				cy.set_bound();
			}*/
		}
		//return;
		cf->setOuterConstraint(constraint_flag, constraint_dir);
		cf->setField();
		cf->initFieldInfo();
		//cf->write_field();
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

		//double s1 = lp.GetV(hl_transfer->to) - v_para;
		//do
		//{
		//	if (visited_f[hl_transfer->left])
		//	{
		//		s1 = lp.GetV(hl_transfer->to) - v_para;
		//		if (s0 * s1 < 0)
		//		{
		//			if (s0 > 0)
		//			{
		//				//poh[0].h = mesh->next_halfedge_handle(hitr.handle());
		//				pohl[0].hl = hl_transfer->next;
		//				pohl[0].c = s0 / (s0 - s1);
		//				distance[0] = s0; distance[1] = s1;
		//				++id[0];
		//			}
		//			else
		//			{
		//				//poh[1].h = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(hitr.handle()));
		//				pohl[1].hl = hl_transfer->next->oppo;
		//				pohl[1].c = s1 / (s1 - s0);
		//				++id[1];
		//			}
		//		}
		//	}
		//	s0 = s1;
		//	hl_transfer = hl_transfer->oppo->next;
		//} while (hl_transfer != hl_begin);

		

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
		//while (h.idx() != poh[1].h.idx())

		//std::ofstream file_writer;
		//file_writer.open("C:\\Users\\123\\Desktop\\eid.txt");
		/*while (hl->id != pohl[1].hl->id)
		{
			int vlid = hl->oppo->prev->from;
			if (!visited_v[vlid])
			{
				return false;
			}
			double s = lp.GetV(vlid) - v_para;
			if (s > 0)
			{
				hl = hl->oppo->next;
				distance[0] = s;
			}
			else
			{
				hl = hl->oppo->prev;
				distance[1] = s;
			}
			planar_loop.emplace_back(hl, distance[0] / (distance[0] - distance[1]));
		}*/
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
				pl.front().c * mesh->point(m4.verticelayers[pl.front().hl->from].v/*mesh->from_vertex_handle(pl.front().h)*/) +
				(1 - pl.front().c) * mesh->point(m4.verticelayers[pl.front().hl->to].v/*mesh->to_vertex_handle(pl.front().h)*/) };
			len[i] = 0;
			len[i] += (pos[0] - pos[1]).norm();
			for (int j = 1; j < pln; ++j)
			{
				pos[(j + 1) & 1] = pl[j].c * mesh->point(m4.verticelayers[pl[j].hl->from].v/*mesh->from_vertex_handle(pl[j].h)*/)
					+ (1 - pl[j].c) * mesh->point(m4.verticelayers[pl[j].hl->to].v/*mesh->to_vertex_handle(pl[j].h)*/);
				len[i] += (pos[0] - pos[1]).norm();
			}
			len[i] += (mesh->point(vertex_set[i]->v) - pl.back().c * mesh->point(m4.verticelayers[pl.back().hl->from].v/*mesh->from_vertex_handle(pl.back().h)*/)
				- (1 - pl.back().c) * mesh->point(m4.verticelayers[pl.back().hl->to].v/*mesh->to_vertex_handle(pl.back().h)*/)).norm();
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
			//double t0 = lp.GetRegularU(mesh->from_vertex_handle(poh.h).idx());
			//double t1 = lp.GetRegularU(mesh->to_vertex_handle(poh.h).idx());
			double t0 = lp.GetRegularU(poh.hl->from);
			double t1 = lp.GetRegularU(poh.hl->to);
			if (fabs(t0 - t1) > 0.5) { if (t0 < t1) { t0 += 1.0; } else { t1 += 1.0; } }
			double u = poh.c * t0 + (1 - poh.c) * t1;
			return u - std::floor(u);
		};
		auto setdata = [&](int i)
		{
			double up[2];
			//up[0] = lp.GetRegularU(minmax_id[i]);
			const auto& pl = lp.all_pl[minmax_id[i]];
			int pln = pl.size();
			up[0] = assembleU(pl.front());
			for (int j = 1; j < pln; ++j)
			{
				up[j & 1] = assembleU(pl[j]);
				if ((up[0] - 0.5) * (up[1] - 0.5) <= 0.0 && fabs(up[0] - up[1]) < 0.5)
				{
					Vec3d p0 = pl[j - 1].c * mesh->point(m4.verticelayers[pl[j - 1].hl->from].v/*mesh->from_vertex_handle(pl[j - 1].h)*/)
						+ (1 - pl[j - 1].c) * mesh->point(m4.verticelayers[pl[j - 1].hl->to].v/*mesh->to_vertex_handle(pl[j - 1].h)*/);
					Vec3d p1 = pl[j].c * mesh->point(m4.verticelayers[pl[j].hl->from].v/*mesh->from_vertex_handle(pl[j].h)*/)
						+ (1 - pl[j].c) * mesh->point(m4.verticelayers[pl[j].hl->to].v/*mesh->to_vertex_handle(pl[j].h)*/);
					double lambda = (0.5 - up[(j + 1) & 1]) / fabs(up[0] - up[1]);
					pos[i] = (1 - lambda) * p0 + lambda * p1;
					//dprint("up", up[0], up[1]);
					return true;
				}
			}
			return false;
		};
		setdata(0);
		setdata(1);
		if ((pos[0] - pos[1]).norm() <= YYSS_FAIRLY_SMALL)
			return YYSS_INFINITE;
		//dprint("setdata:", setdata(0), setdata(1));
		//dprint(minmax_len[0], minmax_len[1]);
		//dprint(pos[0], "\t", pos[1]);
		//static int rrr = 0;
		//++rrr;
		//dprint("rrr", rrr);
		////if (rrr > 18)
		//{
		//	u0point5.push_back(pos[0]);
		//	u0point5.push_back(pos[1]);
		//}
		return (minmax_len[1] - minmax_len[0]) / (pos[0] - pos[1]).norm();
	}

	void LoopGen::AssembleIOVLoopEnergy(info_pair_pq& pq)
	{
		//int vvvid = InfoOnMesh[16389 * 2].energy < InfoOnMesh[16389 * 2 + 1].energy ? 16389 * 2 : 16389 * 2 + 1;
		for (int i = 0; i < mesh->n_vertices(); ++i)
		{
			if (m4.sing_flag[i])
				continue;
			for (int j = 0; j < 2; ++j)
			{
				VertexLayer* vl = &m4.verticelayers[m4.verticemap[i] + j];
				auto& iov = InfoOnMesh[vl->id];
				auto& pl = pls[iov.plid];
				if (pl.empty())
					continue;
				double sum = iov.energy;
				for (auto pohl : pl)
				{
					sum += pohl.c * InfoOnMesh[pohl.hl->from].energy + (1 - pohl.c) * InfoOnMesh[pohl.hl->to].energy;
				}
				//int fromid = pl.front().hl->from;
				//int toid = pl.front().hl->to;
				//InfoOnVertex* iov_transfer = &InfoOnMesh[fromid];
				//double e[2];
				//e[0] = iov_transfer->energy;
				//iov_transfer = &InfoOnMesh[toid];
				//e[1] = iov_transfer->energy;
				//sum += pl.front().c * e[0] + (1 - pl.front().c) * e[1];
				//for (int k = 1; k < pl.size(); ++k)
				//{
				//	int from = pl[k].hl->from;
				//	if (from == fromid)
				//	{
				//		toid = pl[k].hl->to;
				//		iov_transfer = &InfoOnMesh[toid];
				//		e[1] = iov_transfer->energy;
				//	}
				//	else
				//	{
				//		fromid = from;
				//		iov_transfer = &InfoOnMesh[fromid];
				//		e[0] = iov_transfer->energy;
				//	}
				//	sum += pl[k].c * e[0] + (1 - pl[k].c) * e[1];
				//	//sum += iov.pl[j].c * e[0] + (1 - iov.pl[j].c) * e[1];
				//}
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
			/*for (auto vv : mesh->vv_range(vh))
			{
				if (vs_flag[vv.idx()] && !visited[vv.idx()])
				{
					tree.push(vv);
					visited[vv.idx()] = true;
				}
			}*/
		}
		//visited.swap(std::deque<bool>());
		if (count < vertex_set.size())
			return false;
		//check manifold
		BoolVector fs_flag(m4.facelayers.size(), false);
		BoolVector bv_flag(m4.verticelayers.size(), false);
		for (const auto& vl : vertex_set)
		{
			/*for (auto vf : mesh->vf_range(vh))
			{
				if (fs_flag[vf.idx()])
					continue;
				bool flag = true;
				for (auto fv : mesh->fv_range(vf))
				{
					if (!vs_flag[fv.idx()])
					{
						flag = false;
						break;
					}
				}
				if (flag)
					fs_flag[vf.idx()] = true;
			}*/
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
						//ht = ht->prev->oppo;
						ht = ht->next;
					} while (ht != hb);
					if (flag)
						fs_flag[hl_transfer->left] = true;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);

			if (bv_flag[vl->id])
				continue;
			/*for (auto vv : mesh->vv_range(vh))
			{
				if (!vs_flag[vv.idx()])
				{
					bv_flag[vh.idx()] = true;
					break;
				}
			}*/
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
			/*for (auto vv : mesh->vv_range(vh))
			{
				if (!vs_flag[vv.idx()])
					continue;
				HalfedgeHandle he = mesh->find_halfedge(vh, vv);
				if (!fs_flag[mesh->face_handle(he).idx()] && !fs_flag[mesh->face_handle(mesh->opposite_halfedge_handle(he)).idx()])
					return false;
			}*/
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
			//if (!bv_flag[vh.idx()])
			if (!bv_flag[vl->id])
				continue;
			//HalfedgeHandle he = mesh->voh_begin(vh).handle(); 
			auto he = vl->hl;
			/*while (!fs_flag[mesh->face_handle(he).idx()] || fs_flag[mesh->face_handle(mesh->opposite_halfedge_handle(he)).idx()])
			{
				he = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));
			}*/
			while (!fs_flag[he->left] || fs_flag[he->oppo->left])
				he = he->prev->oppo;
			//if (be_flag[he.idx() / 2])
			if (be_flag[he->id / 2])
				continue;
			HalfedgeLayer* hl_transfer = he;
			do
			{
				/*grow_dir[mesh->to_vertex_handle(he_transfer).idx()] = count;
				be_flag[he_transfer.idx() / 2] = true;
				he_transfer = mesh->opposite_halfedge_handle(he_transfer);
				while (!fs_flag[mesh->face_handle(he_transfer).idx()] || fs_flag[mesh->face_handle(mesh->opposite_halfedge_handle(he_transfer)).idx()])
				{
					he_transfer = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he_transfer));
				}*/
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
			//if (grow_dir[vh.idx()] == -1)
				//grow_dir[vh.idx()] = 0;
			if (grow_dir[vl->id] == -1)
				grow_dir[vl->id] = 0;
		}
		return true;
	}

#if 0
	void LoopGen::UpdateIOM()
	{
		for (auto& iov : InfoOnMesh)
		{
			iov.mark.clear();
			iov.pl.swap(PlaneLoop());
		}
		auto& matching = cf->getMatching();
		for (auto eitr = mesh->edges_begin(); eitr != mesh->edges_end(); ++eitr)
		{
			int index = 0;
			auto h = mesh->halfedge_handle(eitr.handle(), 0);
			auto fromvert = mesh->from_vertex_handle(h);
			auto tovert = mesh->to_vertex_handle(h);
			auto ht = mesh->voh_begin(fromvert).handle();
			//逆时针搜索
			while (ht.idx() != h.idx())
			{
				ht = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(ht));
				index += matching[ht.idx()];
			}
			h = mesh->next_halfedge_handle(h);
			ht = mesh->voh_begin(tovert).handle();
			while (ht.idx() != h.idx())
			{
				ht = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(ht));
				index += 4 - matching[ht.idx()];
			}
			index %= 4;

			switch (index)
			{
			case 0:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(2 * tovert.idx(), 0));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(2 * fromvert.idx(), 0));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(2 * tovert.idx() + 1, 0));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(2 * fromvert.idx() + 1, 0));
				break;
			case 1:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(2 * tovert.idx() + 1, 0));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(2 * fromvert.idx() + 1, 0));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(2 * tovert.idx(), 0));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(2 * fromvert.idx(), 0));
				break;
			case 2:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(2 * tovert.idx(), 0));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(2 * fromvert.idx(), 0));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(2 * tovert.idx() + 1, 0));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(2 * fromvert.idx() + 1, 0));
				break;
			case 3:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(2 * tovert.idx() + 1, 0));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(2 * fromvert.idx() + 1, 0));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(2 * tovert.idx(), 0));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(2 * fromvert.idx(), 0));
				break;
			}
		}
		auto& crossfield = cf->getCrossField();
		for (auto& cy : cset.cylinders)
		{
			VertexHandle vh = cy.verticelayers.front();
			HalfedgeHandle heh = mesh->voh_begin(vh);
			FaceHandle fh = mesh->face_handle(heh);
			auto grad = cy.GetUGrad(fh);
			int dir_flag = -1;
			double dir_max = -1.0;
			for (int i = 4 * fh.idx(); i < 4 * fh.idx() + 4; ++i)
			{
				double dot_ = grad[0] * crossfield(0, i) + grad[1] * crossfield(1, i) + grad[2] * crossfield(2, i);
				if (dot_ > dir_max)
				{
					dir_flag = i;
					dir_max = dot_;
				}
			}

			int begin = vh.idx() * 2 + dir_flag % 2;
			auto& ior = cy.info_on_region;
			ior.resize(mesh->n_vertices() * 2, false);
			ior[begin] = true;
			std::queue<InfoOnVertex*> vtree;
			vtree.push(&InfoOnMesh[begin]);
			while (!vtree.empty())
			{
				InfoOnVertex* iov = vtree.front();
				vtree.pop();
				cset.iov_cylinder_map[iov->id].push_back(cy.id);
				for (const auto& iov_nei : iov->mark)
				{
					if (ior[iov_nei.first] || !cy.vertice_flag[iov_nei.first / 2])
						continue;
					vtree.push(&InfoOnMesh[iov_nei.first]);
					ior[iov_nei.first] = true;
				}
			}
		}
	}

	void LoopGen::ProcessOverlap()
	{
		cylinder* base;
		auto merge = [&]()
		{

		};
		auto cut = [&]()
		{

		};
		for (int i = 0; i < cset.cylinders.size(); ++i)
		{
			base = &cset.cylinders[i];
			for (auto tv : base->verticelayers)
			{
				//if(i)
			}
		}
	}

	void LoopGen::LookForCrossingLoop()
	{
		UpdateIOM();
		ProcessOverlap();
		//for (auto& cy : cset.cylinders)
		//{
		//	for (auto he : cy.bounds)
		//	{
		//		auto bv = mesh->to_vertex_handle(he);
		//		//if (cset.vertex_cylinder_map[bv.idx()].size() > 1)
		//			//continue;

		//	}
		//}
	}
#endif

#if 0
	void LoopGen::SetUParaLine(InfoOnVertex& iov, LocalParametrization& lp, std::deque<bool>& visited_v, std::deque<bool>& visited_f)
	{
		iov.pl.clear();
		double u_para = lp.GetRegularU(iov.v.idx());
		if (u_para < 0.1 || u_para > 0.9)
			return;
		iov.pl.clear();
		auto hitr = mesh->voh_begin(iov.v);
		if (!visited_f[hitr.handle().face().idx()])
			return;
		double s0 = lp.GetRegularU(mesh->to_vertex_handle(mesh->next_halfedge_handle(hitr.handle())).idx()) - u_para;
		double distance[2];
		PointOnHalfedge poh[2];
		int id[2] = { 0,0 };
		for (; hitr != mesh->voh_end(iov.v); ++hitr)
		{
			double s1 = lp.GetRegularU(mesh->to_vertex_handle(hitr.handle()).idx()) - u_para;
			if (visited_f[hitr.handle().face().idx()])
			{
				if (s0 * s1 < 0)
				{
					if (s0 > 0)
					{
						poh[0].h = mesh->next_halfedge_handle(hitr.handle());
						poh[0].c = s0 / (s0 - s1);
						distance[0] = s0; distance[1] = s1;
						++id[0];
					}
					else
					{
						poh[1].h = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(hitr.handle()));
						poh[1].c = s1 / (s1 - s0);
						++id[1];
					}
				}
			}
			s0 = s1;
		}
		if (id[0] != 1 || id[1] != 1)
			return;

		auto &planar_loop = iov.pl;
		planar_loop.push_back(poh[0]);
		auto h = poh[0].h;
		while (h.idx() != poh[1].h.idx())
		{
			int vid = mesh->from_vertex_handle(mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(h))).idx();
			/*if (il == 9030)
			{
				dprint(vid);
				dprint("halfedge:", h.idx() / 2);
			}*/
			if (!visited_v[vid])
				return;
			double s = lp.GetRegularU(vid) - u_para;
			if (s > 0)
			{
				h = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h));
				distance[0] = s;
			}
			else
			{
				h = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(h));
				distance[1] = s;
			}
			planar_loop.emplace_back(h, distance[0] / (distance[0] - distance[1]));
		}
	}
#endif

}