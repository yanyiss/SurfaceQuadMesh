#include "LoopGen.h"
#include <omp.h>

#define COMPUTE_NEW_PLANELOOP 0
#define COMPUTE_NEW_ENERGY 0
#define PRINT_DEBUG_INFO 0

namespace LoopGen
{
#pragma region initialization
	void LoopGen::GetPositionFromLoop(const std::vector<VertexHandle>& loop, Eigen::VectorXd xyz[3])
	{
		int n = loop.size() - 1;
		xyz[0].resize(n); xyz[1].resize(n); xyz[2].resize(n);
		for (int i = 0; i < n; ++i)
		{
			auto& pos = mesh->point(loop[i]);
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

	double LoopGen::RefineLoopByPlanarity(std::vector<VertexHandle>& loop, PlaneLoop& planar_loop, int shift)
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
		auto hitr = mesh->voh_begin(loop[0]);
		auto s0 = dis(mesh->point(mesh->to_vertex_handle(mesh->next_halfedge_handle(hitr.handle()))));
		double distance[2];
		PointOnHalfedge poh[2]; int id[2] = { 0, 0 };
		for (; hitr != mesh->voh_end(loop[0]); ++hitr)
		{
			auto s1 = dis(mesh->point(mesh->to_vertex_handle(hitr.handle())));
			if (s0 * s1 < 0)
			{
				/*h[vidmap] = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(hitr.handle()));
				++vidmap;*/
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
		}
		if (id[0] != 1 || id[1] != 1)
		{
			//dprint("error in repairing loop:", loop[0].idx());
			return YYSS_INFINITE;
		}

		//检查出发点到poh[0]的方向与搜索loop的方向是否相同，若不相同，则调换poh[0]和poh[1]
		{
			auto& matching = cf->getMatching();
			//int shift = 0;
			auto h_end = mesh->prev_halfedge_handle(poh[0].h);
			auto h_begin = mesh->voh_begin(mesh->from_vertex_handle(h_end)).handle();
			while (h_begin != h_end)
			{
				h_begin = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h_begin));
				shift += 4 - matching[h_begin.idx()];
			}
			auto& v0 = cf->getCrossField().col(4 * mesh->face_handle(h_end).idx() + (shift % 4));
			auto v1 = poh[0].c * mesh->point(mesh->from_vertex_handle(poh[0].h)) +
				(1 - poh[0].c) * mesh->point(mesh->to_vertex_handle(poh[0].h)) - mesh->point(mesh->from_vertex_handle(h_end));
			if (v0(0) * v1[0] + v0(1) * v1[1] + v0(2) * v1[2] < 0)
			{
				std::swap(poh[0], poh[1]);
				for (int i = 0; i < 4; ++i)
					plane[i] *= -1.0;
				for (int i = 0; i < 2; ++i)
				{
					poh[i].h = mesh->opposite_halfedge_handle(poh[i].h);
					poh[i].c = 1 - poh[i].c;
				}
				distance[0] = dis(mesh->point(mesh->to_vertex_handle(poh[0].h)));
				distance[1] = dis(mesh->point(mesh->from_vertex_handle(poh[0].h)));
			}
		}

		planar_loop.push_back(poh[0]);
		auto h = poh[0].h;
		while (h.idx() != poh[1].h.idx())
		{
			auto s = dis(mesh->point(mesh->from_vertex_handle(mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(h)))));
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
			//dprint(h.idx(),h.idx()/2);
			/*if (h.idx() == 17413)
				break;*/
		}
		return EvaluatePlanarity(xyz, plane);
	}

	bool LoopGen::FieldAligned_PlanarLoop(VertexHandle v, std::vector<VertexHandle>& loop, int shift)
	{
		//从顶点第0个半边左侧面的第0个crossfield的相反方向找路径，最后得到与场方向相同的loop
		int nv = mesh->n_vertices();
		int vid = v.idx();
		std::deque<bool> visited(nv, false);
		std::vector<double> distance(nv, YYSS_INFINITE);
		std::vector<HalfedgeHandle> prev(nv);
		shift %= 2;
		shift += 2;

		struct VertexPQ
		{
			int id;
			int shift;
			double dist;
			int count;
			VertexPQ() {}
			VertexPQ(int id_, int shift_, double dist_, int count_) :id(id_), shift(shift_), dist(dist_), count(count_) {}
			bool operator>(const VertexPQ& x) const { return dist > x.dist; }
		};
		std::priority_queue<VertexPQ, std::vector<VertexPQ>, std::greater<VertexPQ>> pq;

		std::vector<int> count(nv, 0);
		Eigen::Vector3d plane_normal(0, 0, 0);
		auto& crossfield = cf->getCrossField();
		auto& matching = cf->getMatching();
		auto& position = cf->getPosition();

		double halfPI = PI * 0.5;
		double doublePI = PI * 2.0;
		double triple_halfPI = halfPI * 3.0;
		for (auto voh = mesh->voh_begin(v); voh != mesh->voh_end(v); ++voh)
		{
			int fid = mesh->face_handle(voh.handle()).idx();
			plane_normal += crossfield.col(4 * fid + shift + 1);
			if (cf->weight(shift, voh->idx()) < YYSS_INFINITE)
			{
				int toid = mesh->to_vertex_handle(voh.handle()).idx();
				distance[toid] = cf->weight(shift, voh->idx());
				pq.emplace(toid, shift, distance[toid], ++count[toid]);
				prev[toid] = voh.handle();
			}
			shift += matching[voh->idx()]; shift %= 4;
		}
		plane_normal.normalize();

		while (true)
		{
			VertexPQ vert;
			do
			{
				if (pq.empty())
				{ 
					return false;
				}
				vert = pq.top();
				pq.pop();
			} while (vert.count != count[vert.id]);

			//loop.push_back(vert.vidmap); loop.push_back(mesh->from_vertex_handle(prev[vert.vidmap]).idx());
			int fromid = vert.id;
			visited[fromid] = true;
			if (fromid == vid)
			{
				break;
			}
			auto voh = mesh->next_halfedge_handle(prev[fromid]);
			int valence = mesh->valence(mesh->vertex_handle(fromid)) - 1;
			shift = vert.shift;
			for (int i = 0; i < valence; ++i)
			{
				double w = cf->weight(shift, voh.idx());// w *= w;
				if (w < YYSS_INFINITE)
				{
					int toid = mesh->to_vertex_handle(voh).idx();
					//double dot_ = (position.col(toid) - position.col(fromid)).normalized().dot(plane_normal); dot_ *= dot_;
					//w = sqrt(w + dot_);
					//w = sqrt(w * 900 + 1 - w);
					if (distance[fromid] + w < distance[toid])
					{
						distance[toid] = distance[fromid] + w;
						pq.emplace(toid, shift, distance[toid], ++count[toid]);
						prev[toid] = voh;
					}
				}
				shift += matching[voh.idx()]; shift %= 4;
				voh = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(voh));
			}
		}

		loop.clear();
		loop.push_back(v);
		auto prevvert = mesh->from_vertex_handle(prev[vid]);
		while (prevvert.idx() != vid)
		{
			loop.push_back(prevvert);
			prevvert = mesh->from_vertex_handle(prev[prevvert.idx()]);
		}
		loop.push_back(v);
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
		timeRecorder tr;
		InfoOnMesh.resize(mesh->n_vertices() * 2);
		for (int i = 0; i < mesh->n_vertices(); ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				InfoOnMesh[2 * i + j].id = 2 * i + j;
			}
		}
		int nedges = mesh->n_edges();

		tr.tog();
		if (COMPUTE_NEW_PLANELOOP || !ReadPlaneLoop(InfoOnMesh, model_name, mesh))
		{
#pragma omp parallel for
			for (int i = 0; i < mesh->n_vertices(); ++i)
			{
				for (int j = 0; j < 2; ++j)
				{
					std::vector<VertexHandle> loop;
					if (FieldAligned_PlanarLoop(mesh->vertex_handle(i), loop, j))
					{
						RefineLoopByPlanarity(loop, InfoOnMesh[2 * i + j].pl, j);
					}
				}
			}
			WritePlaneLoop(InfoOnMesh, model_name, mesh);
		}
		tr.out("Time of Initializing Planar Loops on All Vertices:");

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
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(2 * tovert.idx() + 1, 1));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(2 * fromvert.idx() + 1, 0));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(2 * tovert.idx(), 0));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(2 * fromvert.idx(), 1));
				break;
			case 2:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(2 * tovert.idx(), 1));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(2 * fromvert.idx(), 1));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(2 * tovert.idx() + 1, 1));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(2 * fromvert.idx() + 1, 1));
				break;
			case 3:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(2 * tovert.idx() + 1, 0));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(2 * fromvert.idx() + 1, 1));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(2 * tovert.idx(), 1));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(2 * fromvert.idx(), 0));
				break;
			}
		}
		dprint("plane loop set done");

		auto assembleLoop = [&](Vec3d &start, PlaneLoop &pl, bool if_forward, Eigen::Matrix3Xd & loop)
		{
			loop.resize(3, 1 + pl.size());
			loop.col(0) << start[0], start[1], start[2];
			int c = 0;
			if (if_forward)
			{
				for (auto itr = pl.begin(); itr != pl.end(); ++itr)
				{
					auto pos = mesh->point(mesh->to_vertex_handle(itr->h)) * (1 - itr->c) + mesh->point(mesh->from_vertex_handle(itr->h)) * itr->c;
					loop.col(++c) << pos[0], pos[1], pos[2];
				}
			}
			else
			{
				for (auto itr = pl.rbegin(); itr != pl.rend(); ++itr)
				{
					auto pos = mesh->point(mesh->to_vertex_handle(itr->h)) * (1 - itr->c) + mesh->point(mesh->from_vertex_handle(itr->h)) * itr->c;
					loop.col(++c) << pos[0], pos[1], pos[2];
				}
			}
		};

		tr.tog();
		similarity_energy.resize(2 * nedges, YYSS_INFINITE);
		if (COMPUTE_NEW_ENERGY || !ReadEnergy(similarity_energy, model_name))
		{
#pragma omp parallel for
			for (int k = 0; k < nedges; ++k)
			{
				auto h = mesh->halfedge_handle(k * 2);
				auto fromvert = mesh->from_vertex_handle(h);
				auto tovert = mesh->to_vertex_handle(h);
				for (int i = 0; i < 2; ++i)
				{
					auto& fl = InfoOnMesh[2 * fromvert.idx() + i];
					//dprint(fl.mark.find(&InfoOnMesh[2 * tovert.idx() + i]) != fl.mark.end());
					auto& tl = fl.mark.find(2 * tovert.idx() + i) != fl.mark.end() ?
						InfoOnMesh[2 * tovert.idx() + i] : InfoOnMesh[2 * tovert.idx() + ((i + 1) % 2)];
					if (fl.pl.empty() || tl.pl.empty())
					{
						similarity_energy[2 * k + i] = YYSS_INFINITE;
						continue;
					}

					int flag = fl.mark[tl.id];
					Eigen::Matrix3Xd loop0, loop1;
					assembleLoop(mesh->point(fromvert), fl.pl, true, loop0);
					assembleLoop(mesh->point(tovert), tl.pl, !flag, loop1);

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
					similarity_energy[2 * k + i] = EvaluateSimilarity(loop0, loop1, u0, id);
					//dprint(k, i, "similarity energy:", e);
				}
			}
			WriteEnergy(similarity_energy, model_name);
		}
		tr.out("time of computing similarity energy:");

		tr.tog();
#pragma omp parallel for
		for (int i = 0; i < InfoOnMesh.size(); ++i)
			InfoOnMesh[i].energy = 0;
#pragma omp parallel for
		for (int i = 0; i < nedges; ++i)
		{
			auto h = mesh->halfedge_handle(i * 2);
			int fromid2 = mesh->from_vertex_handle(h).idx() * 2;
			int toid2 = mesh->to_vertex_handle(h).idx() * 2;
			for (int j = 0; j < 2; ++j)
			{
				auto& fl = InfoOnMesh[fromid2 + j];
				auto& tl = fl.mark.find(toid2 + j) != fl.mark.end() ?
					InfoOnMesh[toid2 + j] : InfoOnMesh[toid2 + ((j + 1) % 2)];
				fl.energy += similarity_energy[i * 2 + j];
				tl.energy += similarity_energy[i * 2 + j];
			}
		}
#pragma omp parallel for
		for (int i = 0; i < InfoOnMesh.size(); ++i)
			InfoOnMesh[i].energy /= mesh->valence(mesh->vertex_handle(i / 2));
		tr.out("time of setting vertex energy:");
	}
#pragma endregion 

	void LoopGen::ConstructInitialRegion(InfoOnVertex* iov, LocalParametrization &lp)
	{
		int nv = mesh->n_vertices();
		int nf = mesh->n_faces();

		auto& pl = iov->pl;
		std::vector<InfoOnVertex*> IOV; 
		IOV.push_back(iov);
		std::vector<std::vector<InfoOnVertex*>> advancing_front[2];
		advancing_front[0].push_back(IOV);
		//advancing_front[1].push_back(IOV);
		std::vector<InfoOnVertex*> hierarchy_vertex[2];
		hierarchy_vertex[0].reserve(pl.size() + 3); //hierarchy_vertex[0].push_back(&iov);
		hierarchy_vertex[1].reserve(pl.size() + 3); //hierarchy_vertex[1].push_back(&iov);
		auto fromid = mesh->from_vertex_handle(pl.front().h).idx();
		auto toid = fromid;
		for (auto vv : iov->mark)
		{
			if (vv.first / 2 == fromid)
			{
				hierarchy_vertex[1].push_back(&InfoOnMesh[vv.first]);
				break;
			}
		}
		//将plane loop上的点加入hierarchy_vertex中
		std::deque<bool> visited_v(mesh->n_vertices(), false);
		visited_v[iov->id / 2] = true;
		visited_v[hierarchy_vertex[1].back()->id / 2] = true;
		for (auto pl_b = pl.begin(); pl_b != pl.end(); ++pl_b)
		{
			auto he = pl_b->h;
			if (mesh->from_vertex_handle(he).idx() == hierarchy_vertex[1].back()->id / 2)
			{
				//hierarchy_vertex[0].push_back(&InfoOnMesh[2*mesh->to_vertex_handle(he).idx()+])
				toid = mesh->to_vertex_handle(he).idx();
				const auto& mark = hierarchy_vertex[1].back()->mark;
				if (!visited_v[toid])
				{
					visited_v[toid] = true;
					toid *= 2;
					hierarchy_vertex[0].push_back(mark.find(toid) != mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
				}
			}
			else
			{
				fromid = mesh->from_vertex_handle(he).idx();
				const auto& mark = hierarchy_vertex[0].back()->mark;
				if (!visited_v[fromid])
				{
					visited_v[fromid] = true;
					fromid *= 2;
					hierarchy_vertex[1].push_back(mark.find(fromid) != mark.end() ? &InfoOnMesh[fromid] : &InfoOnMesh[fromid + 1]);
				}
			}
		}

		//将出发点的1邻域加入hierarchy_vertex中
		if (hierarchy_vertex[0].front() == hierarchy_vertex[0].back())
			hierarchy_vertex[0].pop_back();
		else
		{
			auto h_b = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(mesh->find_halfedge(mesh->vertex_handle(iov->id / 2), mesh->vertex_handle(hierarchy_vertex[0].front()->id / 2))));
			int h_e_idx = mesh->find_halfedge(mesh->vertex_handle(iov->id / 2), mesh->vertex_handle(hierarchy_vertex[0].back()->id / 2)).idx();
			if (h_e_idx != -1)
			{
				while (h_b.idx() != h_e_idx)
				{
					toid = mesh->to_vertex_handle(h_b).idx();
					if (!visited_v[toid])
					{
						visited_v[toid] = true;
						toid *= 2;
						hierarchy_vertex[0].push_back(iov->mark.find(toid) != iov->mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
					}
					//hierarchy_vertex[0].push_back(iov->mark.find(&InfoOnMesh[toid]) != iov->mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
					h_b = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h_b));
				}
			}
		}
		if (hierarchy_vertex[1].front() == hierarchy_vertex[1].back())
			hierarchy_vertex[1].pop_back();
		else
		{
			auto h_b = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(mesh->find_halfedge(mesh->vertex_handle(iov->id / 2), mesh->vertex_handle(hierarchy_vertex[1].back()->id / 2))));
			int h_e_idx = mesh->find_halfedge(mesh->vertex_handle(iov->id / 2), mesh->vertex_handle(hierarchy_vertex[1].front()->id / 2)).idx();
			if (h_e_idx != -1)
			{
				while (h_b.idx() != h_e_idx)
				{
					toid = mesh->to_vertex_handle(h_b).idx();
					if (!visited_v[toid])
					{
						visited_v[toid] = true;
						toid *= 2;
						hierarchy_vertex[1].push_back(iov->mark.find(toid) != iov->mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
					}
					//hierarchy_vertex[1].push_back(iov->mark.find(&InfoOnMesh[toid]) != iov->mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
					h_b = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h_b));
				}
			}
		}

		/*auto h_b = mesh->next_halfedge_handle(mesh->find_halfedge(hierarchy_vertex[0].back()->v, iov->v));
		auto h_e_idx = mesh->find_halfedge(iov->v, hierarchy_vertex[0].front()->v).idx();
		for (; h_b.idx() != h_e_idx; h_b = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h_b)))
		{
			toid = 2 * mesh->to_vertex_handle(h_b).idx();
			hierarchy_vertex[0].push_back(iov->mark.find(&InfoOnMesh[toid]) != iov->mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
		}
		h_b = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(mesh->find_halfedge(iov->v, hierarchy_vertex[1].back()->v)));
		h_e_idx = mesh->find_halfedge(iov->v, hierarchy_vertex[1].front()->v).idx();
		for (; h_b.idx() != h_e_idx; h_b = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h_b)))
		{
			toid = 2 * mesh->to_vertex_handle(h_b).idx();
			hierarchy_vertex[1].push_back(iov->mark.find(&InfoOnMesh[toid]) != iov->mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
		}*/
		//advancing_front[0].push_back(std::move(hierarchy_vertex[0]));
		advancing_front[0][0] = std::move(hierarchy_vertex[0]);
		advancing_front[1].push_back(std::move(hierarchy_vertex[1]));

		//double energy_threshold = 0.2;


		while (true)
		{
			const auto& current_af = advancing_front[0].back();
			std::vector<InfoOnVertex*> hierarchy;
			for (auto caf : current_af)
			{
				for (const auto& iov_nei : caf->mark)
				{
					int inid = iov_nei.first / 2;
					if (visited_v[inid])
						continue;
					if (InfoOnMesh[iov_nei.first].energy > energy_threshold)
						goto target0;
					visited_v[inid] = true;
					hierarchy.push_back(&InfoOnMesh[iov_nei.first]);
				}
			}
			advancing_front[0].push_back(std::move(hierarchy));
		}
	target0:;
		while (true)
		{
			const auto& current_af = advancing_front[1].back();
			std::vector<InfoOnVertex*> hierarchy;
			for (auto caf : current_af)
			{
				for (const auto& iov_nei : caf->mark)
				{
					int inid = iov_nei.first / 2;
					if (visited_v[inid])
						continue;
					if (InfoOnMesh[iov_nei.first].energy > energy_threshold)
						goto target1;
					visited_v[inid] = true;
					hierarchy.push_back(&InfoOnMesh[iov_nei.first]);
				}
			}
			advancing_front[1].push_back(std::move(hierarchy));
		}
	target1:;

		auto& new_vertex = lp.new_vertex;
		auto& new_face = lp.new_face;
		auto& newv_flag = lp.new_v_flag; newv_flag.resize(nv, false);
		auto& newf_flag = lp.new_f_flag; newf_flag.resize(nf, false);
		auto& regionv_flag = lp.region_v_flag;
		auto& regionf_flag = lp.region_f_flag;
		//auto& grow_dir = lp.GetGrowDir();
		for (int i = 0; i < 2; ++i)
		{
			for (auto& ss : advancing_front[i])
			{
				for (auto& tt : ss)
				{
					new_vertex.push_back(mesh->vertex_handle(tt->id / 2));
					//grow_dir[tt->v.idx()] = i;
					newv_flag[tt->id / 2] = true;
				}
			}
		}

		for (auto new_v : new_vertex)
		{
			for (auto vf = mesh->vf_begin(new_v); vf != mesh->vf_end(new_v); ++vf)
			{
				int vf_id = vf->idx();
				if (regionf_flag[vf_id] || newf_flag[vf_id])
					continue;
				for (auto vfv = mesh->fv_begin(vf.handle()); vfv != mesh->fv_end(vf.handle()); ++vfv)
				{
					if (!(regionv_flag[vfv->idx()] || newv_flag[vfv->idx()]))
						goto target2;
				}
				new_face.push_back(vf.handle());
				newf_flag[vf_id] = true;
			target2:;
			}
		}

		std::queue<FaceHandle> face_tree;
		face_tree.push(mesh->face_handle(mesh->voh_begin(mesh->vertex_handle(iov->id / 2)).handle()));
		std::vector<int> ff_id(nf, -1);
		//ff_id[face_tree.front().idx()] = iov == &InfoOnMesh[iov->v.idx() * 2] ? 0 : 1;
		ff_id[face_tree.front().idx()] = iov->id % 2;
		auto& matching = cf->getMatching();
		auto& crossfield = cf->getCrossField();
		std::deque<bool> searched(nf, false);
		searched[face_tree.front().idx()] = true;
		
		auto& x_axis = lp.x_axis; 
		auto& y_axis = lp.y_axis; 
		while (!face_tree.empty())
		{
			auto ft = face_tree.front();
			int ftid = ft.idx();
			x_axis.col(ftid) = crossfield.col(4 * ftid + ff_id[ftid]);
			y_axis.col(ftid) = crossfield.col(4 * ftid + (ff_id[ftid] + 1) % 4);
			face_tree.pop();
			for (auto fh = mesh->fh_begin(ft); fh != mesh->fh_end(ft); ++fh)
			{
				auto oppo_f = mesh->face_handle(mesh->opposite_halfedge_handle(fh)).idx();
				if (searched[oppo_f] || !newf_flag[oppo_f])
					continue;
				searched[oppo_f] = true;
				ff_id[oppo_f] = (ff_id[ftid] + matching[fh->idx()]) % 4;
				face_tree.push(mesh->face_handle(oppo_f));
			}
		}
	}

	void LoopGen::AssembleSimilarityAngle(VertexHandle v, Eigen::VectorXd& sa, LocalParametrization& lp, int loop_fragment_num)
	{
		auto setData = [&](double t0, double t1, double c, Eigen::Matrix3Xd& fragment, double fromu, double& tou, OpenMesh::Vec3d& frompos, OpenMesh::Vec3d& topos)
		{
			//dprint();
			//dprint("t", t0, t1);
			if (fabs(t0 - t1) > 0.5) { if (t0 < t1) { t0 += 1.0; } else { t1 += 1.0; } }
			tou = c * t0 + (1 - c) * t1;
			tou -= std::floor(tou);
			//dprint("ftu", fromu, tou);
			if (fromu > tou && fabs(fromu - tou) < 0.5) return;
			int u0 = std::floor(fromu * loop_fragment_num);
			int u1 = std::floor(tou * loop_fragment_num);
			//dprint("u", u0, u1);
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
		frompos = mesh->point(v);
		fromu = lp.GetRegularU(v.idx());
		//for (auto& pl : v_ptr->pl)
		for (auto& pl : lp.all_pl[v.idx()])
		{
			topos = pl.c * mesh->point(mesh->from_vertex_handle(pl.h)) + (1 - pl.c) * mesh->point(mesh->to_vertex_handle(pl.h));
			setData(lp.GetU(mesh->from_vertex_handle(pl.h).idx()), lp.GetU(mesh->to_vertex_handle(pl.h).idx()), pl.c, fragment, fromu, tou, frompos, topos);
			frompos = topos;
			fromu = tou;
		}
		topos = mesh->point(v);
		tou = lp.GetRegularU(v.idx());
		setData(tou, tou, 0, fragment, fromu, tou, frompos, topos);
		/*for (int i = 0; i < fragment.cols(); ++i)
		{
			dprint(fragment(0, i), fragment(1, i), fragment(2, i));
		}*/
		sa.resize(loop_fragment_num * (loop_fragment_num + 1) / 2); sa.setZero();
		int count = 0;
		for (int i = 0; i < loop_fragment_num; ++i)
		{
			for (int j = i + 1; j < loop_fragment_num; ++j)
			{
				sa(count++) = conservativeArcCos(fragment.col(i).dot(fragment.col(j)));
			}
		}
		/*for (int i = 0; i < sa.size(); ++i)
			if (fabs(sa(i)) < YYSS_FAIRLY_SMALL)
				sa(i) = YYSS_FAIRLY_SMALL;*/
	}

	bool LoopGen::SpreadSubRegion(LocalParametrization& lp, bool grow_flag[2])
	{
		auto& newv_flag = lp.new_v_flag;
		auto& newf_flag = lp.new_f_flag;
		auto& new_vertex = lp.new_vertex;
		auto& new_face = lp.new_face;
		std::deque<bool> visited_v = lp.region_v_flag;
		std::deque<bool> visited_f = lp.region_f_flag;

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

		//提取含有完整loop的区域
		if (lp.region_vertex.size() == 1)
			RefineLoopByParametrization(lp.region_vertex.front(), lp, visited_v, visited_f);
		std::vector<VertexHandle> vertex_cache;
		std::deque<bool> vertex_cache_flag = lp.region_v_flag;
		for (auto new_v : new_vertex)
		{
			////////////////////////注意！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
			//这里取得InfoOnVertex是随意取的，后续需要更改
			if (RefineLoopByParametrization(new_v, lp, visited_v, visited_f))
			{
				vertex_cache_flag[new_v.idx()] = true;
				vertex_cache.push_back(new_v);
				//SetUParaLine(InfoOnMesh[2 * newid + 1], lp, visited_v, visited_f);
			}
		}
#if PRINT_DEBUG_INFO
		dprint("提取含有完整loop的区域");
#endif
		if (lp.region_face.empty())
		{
			for (auto vv : mesh->vv_range(lp.region_vertex.front()))
			{
				if (!vertex_cache_flag[vv.idx()])
					return false;
			}
			if (!CheckTopology(vertex_cache, vertex_cache_flag, lp.grow_dir))
				return false;
		}
		if (vertex_cache.empty())
			return false;

		auto& regionv_flag = lp.region_v_flag;
		auto& regionf_flag = lp.region_f_flag;
		auto& grow_dir = lp.grow_dir;

		//检测相似性能量
		auto& normal_similarity_angle = lp.normal_similarity_angle;
		//int loop_fragment_num = InfoOnMesh[2 * lp.GetRegionVertex().front().idx()].pl.size();
		auto& all_pl = lp.all_pl;
		int loop_fragment_num = all_pl[lp.region_vertex.front().idx()].size();
		if (vertex_cache_flag[lp.region_vertex.front().idx()] && !lp.has_nsa)
		{
			lp.has_nsa = true;
			AssembleSimilarityAngle(lp.region_vertex.front(), normal_similarity_angle, lp, loop_fragment_num);
		}

		std::deque<bool> if_similarity_energy_low;
		if (lp.region_vertex.size() > 1)
		{
			if_similarity_energy_low.resize(mesh->n_vertices(), false);
			int exceed[2] = { 0,0 };

			static omp_lock_t lock;
			omp_init_lock(&lock);
#pragma omp parallel for
			for (int i = 0; i < vertex_cache.size(); ++i)
			{
				int new_id = vertex_cache[i].idx();
				int grow_id = grow_dir[new_id];
				if (!grow_flag[grow_id])
					continue;
				Eigen::VectorXd similarity_angle;
				AssembleSimilarityAngle(vertex_cache[i], similarity_angle, lp, loop_fragment_num);
				double sum = 0;
				for (int j = 0; j < similarity_angle.size(); ++j) {
					//dprint(j);
					//dprint(similarity_angle(j), normal_similarity_angle(j), inverse_normal_similarity_angle(j));
					//sum += std::min(100.0, fabs(normal_similarity_angle(j) / similarity_angle(j)) + fabs(similarity_angle(j) * inverse_normal_similarity_angle(j)) - 2.0);
					sum += fabs(normal_similarity_angle(j) - similarity_angle(j));
				}
				//dprint(sum / similarity_angle.size());
				if (sum < energy_threshold * similarity_angle.size())
				{
					if_similarity_energy_low[new_id] = true;
					//omp_set_lock(&lock);
					//if (exceed[grow_id] > 2)
					//{
					//	//if_similarity_energy_low[new_id] = false;
					//	grow_flag[grow_id] = false;
					//}
					//++exceed[grow_id];
					//omp_unset_lock(&lock);
				}
			}
			omp_destroy_lock(&lock);
			//normal_similarity_angle = ( normal_similarity_angle * lp.GetRegionVertex().size() + similarity_angle_sum) / (lp.GetRegionVertex().size() + vertex_cache.size());

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
				{
					++exceed[grow_dir[vc.idx()]];
				}
				if (exceed[grow_dir[vc.idx()]] > 2)
				{
					grow_flag[grow_dir[vc.idx()]] = false;
				}
			}
	    }
		else
			if_similarity_energy_low.resize(mesh->n_vertices(), true);
		//dprint("grow_flag", grow_flag[0], grow_flag[1]);
		//grow_flag[1] = false;
		for (int i = 0; i < 2; ++i)
		{
			if (grow_flag[i])
			{
				double grad = LoopLenGrad(vertex_cache, lp, vertex_cache_flag, i == 0 ? false : true);
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
			//if (grow_flag[grow_dir[new_v.idx()]])
			if (if_similarity_energy_low[new_v.idx()])
			{
				region_vertex.push_back(new_v);
				int newid = new_v.idx();
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
			for (auto fv = mesh->fv_begin(new_f); fv != mesh->fv_end(new_f); ++fv)
			{
				if (!regionv_flag[fv->idx()])
					goto target1;
			}
			int newid = new_f.idx();
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
		//int extend_layer = 3;
		newv_flag.resize(mesh->n_vertices(), false);
		new_vertex.clear();
		for (int i = begin_; i < count; ++i)
		{
			auto rvi = region_vertex[i];
			int growid = grow_dir[rvi.idx()];
			if (!grow_flag[growid])
				continue;
			for (auto vv = mesh->vv_begin(rvi); vv != mesh->vv_end(rvi); ++vv)
			{
				int vvid = vv->idx();
				if (regionv_flag[vvid] || newv_flag[vvid])
					continue;
				newv_flag[vvid] = true;
				new_vertex.push_back(vv.handle());
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
				//dprint(growid);
				for (auto vv = mesh->vv_begin(rvj); vv != mesh->vv_end(rvj); ++vv)
				{
					int vvid = vv->idx();
					if (regionv_flag[vvid] || newv_flag[vvid])
						continue;
					newv_flag[vvid] = true;
					new_vertex.push_back(vv.handle());
					grow_dir[vvid] = growid;
				}
			}
			begin_ = end_;
		}

		newf_flag.resize(mesh->n_faces(), false);
		new_face.clear();
		for (auto new_v : new_vertex)
		{
			int newid = new_v.idx();
			for (auto vf = mesh->vf_begin(new_v); vf != mesh->vf_end(new_v); ++vf)
			{
				int vfid = vf->idx();
				if (regionf_flag[vfid] || newf_flag[vfid])
					continue;
				for (auto vfv = mesh->fv_begin(vf.handle()); vfv != mesh->fv_end(vf.handle()); ++vfv)
				{
					if (!regionv_flag[vfv->idx()] && !newv_flag[vfv->idx()])
						goto target2;
				}
				newf_flag[vfid] = true;
				new_face.push_back(vf.handle());
			target2:;
			}
		}
#if PRINT_DEBUG_INFO
		dprint("扩展advancing_front");
#endif
		//sub_vertex = new_vertex;
		//sub_face = new_face;
		
		//优化新区域的场
		ResetLocalField(lp, new_face, newf_flag, regionf_flag);

#if PRINT_DEBUG_INFO
		dprint("优化新区域的场");
#endif
		return true;
	}

	void LoopGen::ResetLocalField(LocalParametrization& lp, std::vector<FaceHandle>& opt_face, std::deque<bool>& opt_flag, std::deque<bool>& constraint_flag)
	{
		typedef std::complex<double> COMPLEX;
		std::vector<Eigen::Triplet<COMPLEX>> triple;
		int count = 0;
		std::vector<int> idmap(mesh->n_faces());
		for (auto f : opt_face)
		{
			idmap[f.idx()] = count++;
		}
		Eigen::VectorXcd b(3 * opt_face.size()); b.setZero();
		count = 0;
		auto& position = cf->getPosition();
		auto& faceBase = cf->getFaceBase();
		auto& x_axis = lp.x_axis;
		auto& y_axis = lp.y_axis;
		for (auto f : opt_face)
		{
			int fid = f.idx();
			for (auto fh = mesh->fh_begin(f); fh != mesh->fh_end(f); ++fh)
			{
				if (fh.handle().edge().is_boundary())
					continue;
				auto gid = fh.handle().opp().face().idx();
				if (!(opt_flag[gid] || constraint_flag[gid]) || (opt_flag[gid] && fid < gid))
					continue;
				auto ev = (position.col(fh.handle().to().idx()) - position.col(fh.handle().from().idx())).normalized();
				COMPLEX e_f = COMPLEX(ev.dot(faceBase.col(fid * 2)), -ev.dot(faceBase.col(fid * 2 + 1)));
				COMPLEX e_g = COMPLEX(ev.dot(faceBase.col(gid * 2)), -ev.dot(faceBase.col(gid * 2 + 1)));
				if (opt_flag[fid])
				{
					triple.emplace_back(count, idmap[fid], e_f);
				}
				else
				{
					COMPLEX dir = COMPLEX(x_axis.col(fid).dot(faceBase.col(2 * fid)), x_axis.col(fid).dot(faceBase.col(2 * fid + 1)));
					b(count) -= e_f * dir;
				}
				if (opt_flag[gid])
				{
					triple.emplace_back(count, idmap[gid], -e_g);
				}
				else
				{
					COMPLEX dir = COMPLEX(x_axis.col(gid).dot(faceBase.col(2 * gid)), x_axis.col(gid).dot(faceBase.col(2 * gid + 1)));
					b(count) = +e_g * dir;
				}
				++count;
			}
		}
		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<COMPLEX>> solver;
		Eigen::SparseMatrix<COMPLEX> A(count, opt_face.size());
		A.setFromTriplets(triple.begin(), triple.end());
		solver.compute(A);
		//dprint(solver.info());
		b = solver.solve(b.head(count));
		/*Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> slu;
		Eigen::SparseMatrix<COMPLEX> A(count, opt_face.size());
		A.setFromTriplets(triple.begin(), triple.end());
		Eigen::SparseMatrix<COMPLEX> AT = A.adjoint();
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> solver;
		solver.compute(AT * A);
		b = solver.solve(AT * b.head(count));*/
		//auto& crossfield = cf->getCrossField();
		for (auto f : opt_face)
		{
			int fid = f.idx();
			double theta = std::arg(b(idmap[fid]));
			x_axis.col(fid) = faceBase.col(fid * 2) * cos(theta) + faceBase.col(fid * 2 + 1) * sin(theta);
			y_axis.col(fid) = faceBase.col(fid * 2) * cos(theta + 0.5 * PI) + faceBase.col(fid * 2 + 1) * sin(theta + 0.5 * PI);
		}
	}

	bool LoopGen::ConstructRegionCut(InfoOnVertex* iov, std::deque<bool>& visited, std::vector<VertexHandle>& cut)
	{
		cut.clear();
		cut.push_back(mesh->vertex_handle(iov->id / 2));
		const auto& matching = cf->getMatching();
		//int shift = iov != &InfoOnMesh[iov->v.idx() * 2] ? 2 : 1;
		int shift = (iov->id % 2) + 1;
		int itertimes = 0;
		for (int s = shift; itertimes < 2; shift += 2, s = shift % 4, ++itertimes)
		{
			HalfedgeHandle prevhe; prevhe.invalidate();
			int smark = s;
			HalfedgeHandle hb = mesh->voh_begin(mesh->vertex_handle(iov->id / 2));
			int hb_idx = hb.idx();
			while (true)
			{
				double w = YYSS_INFINITE;
				do
				{
					if (cf->weight(s, hb.idx()) < w)
					{
						w = cf->weight(s, hb.idx());
						prevhe = hb;
						smark = s;
					}
					s += matching[hb.idx()]; s %= 4;
					hb = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(hb));
				} while (hb.idx() != hb_idx);
				if (!prevhe.is_valid())
				{
					return false;
				}
				if (!visited[mesh->to_vertex_handle(prevhe).idx()])
					break;
				cut.push_back(mesh->to_vertex_handle(prevhe));
				hb_idx = mesh->opposite_halfedge_handle(prevhe).idx();
				s = smark;
				hb = mesh->next_halfedge_handle(prevhe);
			}
			std::reverse(cut.begin(), cut.end());
		}
#if PRINT_DEBUG_INFO
		dprint("计算cut");
#endif
		return true;
	}

	void LoopGen::OptimizeLoop()
	{
		timeRecorder tr;
		info_pair_pq pq;
		AssembleIOVLoopEnergy(pq);
		//info_pair_pq fe; fe.emplace(InfoOnMesh[9051 * 2].energy < InfoOnMesh[9051 * 2 + 1].energy ? &InfoOnMesh[9051 * 2] : &InfoOnMesh[9051 * 2 + 1], 0);
		//pq = fe;

		std::deque<bool> ifset_flag(2 * mesh->n_vertices(), false);
		std::deque<bool> constraint_flag(mesh->n_faces(), false);
		Eigen::Matrix3Xd constraint_dir(3, mesh->n_faces());
		all_plane_loop.resize(mesh->n_vertices());
		optimized_face_flag.resize(mesh->n_faces(), false);
		optimized_vert_flag.resize(mesh->n_vertices(), false);
		bound_edge_flag.resize(mesh->n_edges(), false);
		//cset.vertex_cylinder_map.resize(mesh->n_vertices());
		while (true)
		{
			info_pair ip;
			ip = pq.top(); pq.pop();
			int vid = ip.iov->id / 2;
			while (ifset_flag[ip.iov == &InfoOnMesh[vid * 2] ? vid * 2 : vid * 2 + 1])
			{
				if (pq.empty())
				{
					ip.energy = YYSS_INFINITE;
					break;
				}
				ip = pq.top(); pq.pop();
				vid = ip.iov->id / 2;
			}
			if (ip.energy > energy_threshold)
				break;
			//if (vid == 20683)
				//break;
			//dprint("seed vertex:", seed_vertex.size() - 1, vid, ip.energy);
			/*dprint("\n\n\n");
			dprint("seed vertex:", vid);
			dprint("similarigy energy:", ip.energy);*/
			LocalParametrization lp(*mesh, mesh->vertex_handle(ip.iov->id / 2));
			ConstructInitialRegion(ip.iov, lp);
			bool grow_flag[2] = { true, true };
			do
			{
				std::deque<bool> visited_v = lp.region_v_flag;
				for (auto& ver : lp.new_vertex)
					visited_v[ver.idx()] = true;
				if (!ConstructRegionCut(ip.iov, visited_v, lp.cut))
					break;
				lp.run(cf->getNormal());
			} while (SpreadSubRegion(lp, grow_flag));

			if (lp.region_vertex.size() > 1)
			{
				seed_vertex.push_back(mesh->vertex_handle(ip.iov->id / 2));
				dprint("seed vertex:", seed_vertex.size() - 1, vid, ip.energy);
			}
			auto& rvf = lp.region_v_flag;
			std::queue<InfoOnVertex*> vtree;
			vtree.push(ip.iov);
			std::deque<bool> if_visited(mesh->n_vertices(), false);
			if_visited[ip.iov->id / 2] = true;
			while (!vtree.empty())
			{
				InfoOnVertex* c = vtree.front();
				vtree.pop();
				//int cid = c->v.idx() * 2;
				//ifset_flag[c == &InfoOnMesh[cid] ? cid : cid + 1] = true;
				ifset_flag[c->id] = true;
				for (const auto& iov_nei : c->mark)
				{
					if (rvf[iov_nei.first / 2] && !if_visited[iov_nei.first / 2])
					{
						vtree.push(&InfoOnMesh[iov_nei.first]);
						if_visited[iov_nei.first / 2] = true;
					}
				}
			}
			{
				int count = 0;
				for (auto isf : ifset_flag)
					if (isf)
						++count;
				//dprint("count:", count);
			}

			for (auto rf : lp.region_face)
			{
				constraint_flag[rf.idx()] = true;
				constraint_dir.col(rf.idx()) = lp.x_axis.col(rf.idx());
			}

			//all_plane_loop.insert(all_plane_loop.end(), lp.GetAllPL().begin(), lp.GetAllPL().end());
			//region_vertex.insert(region_vertex.end(), lp.GetRegionVertex().begin(), lp.GetRegionVertex().end());
			auto& region_face = lp.region_face;
			auto& apl = lp.all_pl;
			for (int i = 0; i < apl.size(); ++i)
			{
				if (apl[i].empty())
					continue;
				all_plane_loop[i] = std::move(apl[i]);
			}
			for (auto& rf : lp.region_face)
				optimized_face_flag[rf.idx()] = true;
			if (lp.region_vertex.size() > 1)
				for (auto& rv : lp.region_vertex)
					optimized_vert_flag[rv.idx()] = true;
			for (auto ff : lp.region_face)
				for (auto fh : mesh->fh_range(ff))
					if (!lp.region_f_flag[mesh->face_handle(mesh->opposite_halfedge_handle(fh)).idx()])
						bound_edge_flag[fh.idx() / 2] = true;

			//cylinder cy;
			//cy.mesh = mesh;
			//cy.id = cset.cylinders.size();
			//cy.vertices = std::move(lp.region_vertex);
			//cy.faces = std::move(lp.region_face);
			//cy.vertice_flag = std::move(lp.region_v_flag);
			//cy.face_flag = std::move(lp.region_f_flag);
			//cy.cut_v_flag = std::move(lp.cutv_flag);
			//cy.cut_f_flag = std::move(lp.cutf_flag);
			//cy.vidmap = std::move(lp.vidmap);
			//cy.uv[0] = std::move(lp.GetU());
			//cy.uv[1] = std::move(lp.GetV());
			////for (auto cyv : cy.vertices)
			//	//cset.vertex_cylinder_map[cyv.idx()].push_back(cy.id);
			//cy.SetBound();
			//cset.cylinders.push_back(std::move(cy));
			//cset.push_back(cy);
			//break;
		}
		cf->setOuterConstraint(constraint_flag, constraint_dir);
		cf->setField();
		cf->initFieldInfo();
		tr.out("repair field:");
		dprint("complete");
	}

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
			VertexHandle vh = cy.vertices.front();
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
			for (auto tv : base->vertices)
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

	bool LoopGen::RefineLoopByParametrization(VertexHandle v, LocalParametrization& lp, std::deque<bool>& visited_v, std::deque<bool>& visited_f)
	{
		double v_para = lp.GetV(v.idx());
		auto hitr = mesh->voh_begin(v);
		if (!visited_f[hitr.handle().face().idx()])
			return false;
		double s0 = lp.GetV(mesh->to_vertex_handle(mesh->next_halfedge_handle(hitr.handle())).idx()) - v_para;
		double distance[2];
		PointOnHalfedge poh[2];
		int id[2] = { 0,0 };
		for (; hitr != mesh->voh_end(v); ++hitr)
		{
			/*if (il == 9030)
			{
				dprint(hitr->idx() / 2);
			}*/
			double s1 = lp.GetV(mesh->to_vertex_handle(hitr.handle()).idx()) - v_para;
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
			return false;

		PlaneLoop planar_loop;
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
				return false;
			double s = lp.GetV(vid) - v_para;
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
		lp.all_pl[v.idx()] = std::move(planar_loop);
		//iov.pl = std::move(planar_loop);
		return true;
	}

	double LoopGen::LoopLenGrad(std::vector<VertexHandle>& vertex_set, LocalParametrization& lp, std::deque<bool>& vertex_flag, int growDir)
	{
		int vsn = vertex_set.size();
		std::vector<double> len(vsn, -1.0);
		auto& grow_dir = lp.grow_dir;
#pragma omp parallel for
		for (int i = 0; i < vsn; ++i)
		{
			if (grow_dir[vertex_set[i].idx()] != growDir)
				continue;
			/*bool if_compute = false;
			for (auto vv : mesh->vv_range(vertex_set[i]))
			{
				if (!vertex_flag[vertex_set[i].idx()])
					if_compute = true;
			}
			if (!if_compute)
				continue;*/
			const auto& pl = lp.all_pl[vertex_set[i].idx()];
			int pln = pl.size();
			Vec3d pos[2] = { mesh->point(vertex_set[i]), 
				pl.front().c * mesh->point(mesh->from_vertex_handle(pl.front().h)) + (1 - pl.front().c) * mesh->point(mesh->to_vertex_handle(pl.front().h)) };
			len[i] = 0;
			len[i] += (pos[0] - pos[1]).norm();
			for (int j = 1; j < pln; ++j)
			{
				pos[(j + 1) & 1] = pl[j].c * mesh->point(mesh->from_vertex_handle(pl[j].h)) + (1 - pl[j].c) * mesh->point(mesh->to_vertex_handle(pl[j].h));
				len[i] += (pos[0] - pos[1]).norm();
			}
			len[i] += (mesh->point(vertex_set[i]) - pl.back().c * mesh->point(mesh->from_vertex_handle(pl.back().h))
				- (1 - pl.back().c) * mesh->point(mesh->to_vertex_handle(pl.back().h))).norm();
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
				minmax_id[0] = vertex_set[i].idx();
			}
			if (len[i] > minmax_len[1])
			{
				minmax_len[1] = len[i];
				minmax_id[1] = vertex_set[i].idx();
			}
		}
		if (minmax_id[0] < 0 || minmax_id[1] < 0)
			return YYSS_INFINITE;
		if (minmax_id[0] == minmax_id[1])
			return YYSS_INFINITE;
		//计算 minmax_id 中两条 loop 在 u = 0.5 处的点的距离
		Vec3d pos[2];
		auto assembleU = [&](const PointOnHalfedge& poh)
		{
			double t0 = lp.GetRegularU(mesh->from_vertex_handle(poh.h).idx());
			double t1 = lp.GetRegularU(mesh->to_vertex_handle(poh.h).idx());
			//dprint("t0t1", t0, t1);
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
					Vec3d p0 = pl[j - 1].c * mesh->point(mesh->from_vertex_handle(pl[j - 1].h))
						+ (1 - pl[j - 1].c) * mesh->point(mesh->to_vertex_handle(pl[j - 1].h));
					Vec3d p1 = pl[j].c * mesh->point(mesh->from_vertex_handle(pl[j].h))
						+ (1 - pl[j].c) * mesh->point(mesh->to_vertex_handle(pl[j].h));
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

		static omp_lock_t lock;
		omp_init_lock(&lock);
#pragma omp parallel for
		for (int i = 0; i < InfoOnMesh.size(); ++i)
		//int i = vvvid;
		{
			//dprint(i);
			auto& iov = InfoOnMesh[i];
			if (iov.pl.empty())
			{
				omp_set_lock(&lock);
				pq.emplace(&iov, YYSS_INFINITE);
				omp_unset_lock(&lock);
				continue;
			}
			double sum = iov.energy;
			//int vid = mesh->from_vertex_handle(iov.pl.front().h).idx();
			//e[0] = &iov == &InfoOnMesh[2 * vid] ? iov.energy : InfoOnMesh[2 * vid + 1].energy;
			int fromid = mesh->from_vertex_handle(iov.pl.front().h).idx();
			int toid = mesh->to_vertex_handle(iov.pl.front().h).idx();
			InfoOnVertex* iov_transfer = &InfoOnMesh[iov.mark.find(2 * fromid) != iov.mark.end() ? 2 * fromid : 2 * fromid + 1];
			double e[2];
			e[0] = iov_transfer->energy;
			iov_transfer = &InfoOnMesh[iov.mark.find(2 * toid) != iov.mark.end() ? 2 * toid : 2 * toid + 1];
			e[1] = iov_transfer->energy;
			sum += iov.pl.front().c * e[0] + (1 - iov.pl.front().c) * e[1];
			for (int j = 1; j < iov.pl.size(); ++j)
			{
				int from = mesh->from_vertex_handle(iov.pl[j].h).idx();
				if (from == fromid)
				{
					toid = mesh->to_vertex_handle(iov.pl[j].h).idx();
					iov_transfer = &InfoOnMesh[iov_transfer->mark.find(2 * toid) != iov_transfer->mark.end() ? 2 * toid : 2 * toid + 1];
					e[1] = iov_transfer->energy;
				}
				else
				{
					fromid = from;
					iov_transfer = &InfoOnMesh[iov_transfer->mark.find(2 * fromid) != iov_transfer->mark.end() ? 2 * fromid : 2 * fromid + 1];
					e[0] = iov_transfer->energy;
				}
				sum += iov.pl[j].c * e[0] + (1 - iov.pl[j].c) * e[1];
				//dprint(e[0], e[1], iov.pl[j].c * e[0] + (1 - iov.pl[j].c) * e[1], sum);
			}

			omp_set_lock(&lock);
			pq.emplace(&iov, sum / (1 + iov.pl.size()));
			omp_unset_lock(&lock);
			//dprint(sum / (1 + iov.pl.size()));
		}
		omp_destroy_lock(&lock);
#if PRINT_DEBUG_INFO
		dprint("计算圈上的平均能量");
#endif
	}

	bool LoopGen::CheckTopology(std::vector<VertexHandle>& vertex_set, std::deque<bool>& vs_flag, std::vector<int>& grow_dir)
	{
		//check if empty
		if (vertex_set.empty())
			return false;
		//check connectivity
		std::deque<bool> visited(mesh->n_vertices(), false);
		std::queue<VertexHandle> tree;
		tree.push(vertex_set.front());
		visited[vertex_set.front().idx()] = true;
		int count = 0;
		while (!tree.empty())
		{
			auto vh = tree.front();
			tree.pop();
			++count;
			for (auto vv : mesh->vv_range(vh))
			{
				if (vs_flag[vv.idx()] && !visited[vv.idx()])
				{
					tree.push(vv);
					visited[vv.idx()] = true;
				}
			}
		}
		visited.swap(std::deque<bool>());
		if (count < vertex_set.size())
			return false;
		//check manifold
		std::deque<bool> fs_flag(mesh->n_faces(), false);
		std::deque<bool> bv_flag(mesh->n_vertices(), false);
		for (const auto& vh : vertex_set)
		{
			for (auto vf : mesh->vf_range(vh))
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
			}
			if (bv_flag[vh.idx()])
				continue;
			for (auto vv : mesh->vv_range(vh))
			{
				if (!vs_flag[vv.idx()])
				{
					bv_flag[vh.idx()] = true;
					break;
				}
			}
		}
		for (const auto& vh : vertex_set)
		{
			for (auto vv : mesh->vv_range(vh))
			{
				if (!vs_flag[vv.idx()])
					continue;
				HalfedgeHandle he = mesh->find_halfedge(vh, vv);
				if (!fs_flag[mesh->face_handle(he).idx()] && !fs_flag[mesh->face_handle(mesh->opposite_halfedge_handle(he)).idx()])
					return false;
			}
		}
		//check boundary number and set grow_dir in passing
		std::deque<bool> be_flag(mesh->n_edges(), false);
		count = 0;
		for (const auto& vh : vertex_set)
		{
			if (!bv_flag[vh.idx()])
				continue;
			HalfedgeHandle he = mesh->voh_begin(vh).handle(); 
			while (!fs_flag[mesh->face_handle(he).idx()] || fs_flag[mesh->face_handle(mesh->opposite_halfedge_handle(he)).idx()])
			{
				he = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));
			}
			if (be_flag[he.idx() / 2])
				continue;
			HalfedgeHandle he_transfer = he;
			do
			{
				grow_dir[mesh->to_vertex_handle(he_transfer).idx()] = count;
				be_flag[he_transfer.idx() / 2] = true;
				he_transfer = mesh->opposite_halfedge_handle(he_transfer);
				while (!fs_flag[mesh->face_handle(he_transfer).idx()] || fs_flag[mesh->face_handle(mesh->opposite_halfedge_handle(he_transfer)).idx()])
				{
					he_transfer = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he_transfer));
				}
			} while (he != he_transfer);
			++count;
			if (count > 2)
				return false;
		}
		if (count < 2)
			return false;
		for (const auto& vh : vertex_set)
		{
			if (grow_dir[vh.idx()] == -1)
				grow_dir[vh.idx()] = 0;
		}
		return true;
	}


















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