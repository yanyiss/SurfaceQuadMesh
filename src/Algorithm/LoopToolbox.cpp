#include "LoopToolbox.h"
#include <omp.h>

#define YYSS_INFINITE 1.0e12
#define YYSS_FAIRLY_SMALL 1.0e-3
#define COMPUTE_NEW_PLANELOOP 0
#define COMPUTE_NEW_ENERGY 0
#define PRINT_DEBUG_INFO 1

namespace LoopGen
{
	void LocalParametrization::run()
	{
		//const auto& matching = cf->getMatching();
		int nf = mesh->n_faces();

		//标记与cut相关的顶点和面，计算装配矩阵需要的数据
		int nv = mesh->n_vertices();
		cutv_flag.resize(nv, false);
		cutf_flag.resize(nf, false);
		std::vector<std::map<VertexHandle, std::pair<bool, Eigen::Vector3d>>> info(nf);
		{
			for (auto c : cut)
				cutv_flag[c.idx()] = true;

			const auto& normal = cf->getNormal();
			HalfedgeHandle he = mesh->find_halfedge(cut[0], cut[1]);
			FaceHandle f = mesh->face_handle(he);
			int fid = f.idx();
			auto calc_vector = [&](bool flag)
			{
				double inv_area = 1.0 / (2 * mesh->calc_face_area(f));
				auto vi = mesh->to_vertex_handle(he);
				auto ev = mesh->calc_edge_vector(mesh->prev_halfedge_handle(he));
				info[fid].insert(std::make_pair(vi, std::make_pair(flag && cutv_flag[vi.idx()], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
				vi = mesh->to_vertex_handle(mesh->next_halfedge_handle(he));
				ev = mesh->calc_edge_vector(he);
				info[fid].insert(std::make_pair(vi, std::make_pair(flag && cutv_flag[vi.idx()], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
				vi = mesh->from_vertex_handle(he);
				ev = mesh->calc_edge_vector(mesh->next_halfedge_handle(he));
				info[fid].insert(std::make_pair(vi, std::make_pair(flag && cutv_flag[vi.idx()], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
			};

			while (new_f_flag[fid])
			{
				cutf_flag[fid] = true;
				calc_vector(true);
				he = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));
				f = mesh->face_handle(he);
				fid = f.idx();
			}
			for (int i = 1; i < cut.size() - 1; ++i)
			{
				he = mesh->find_halfedge(cut[i], cut[i + 1]);
				f = mesh->face_handle(he);
				fid = f.idx();
				do
				{
					if (info[fid].empty() && new_f_flag[fid])
					{
						cutf_flag[fid] = true;
						calc_vector(true);
					}
					he = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));
					f = mesh->face_handle(he);
					fid = f.idx();
				} while (mesh->to_vertex_handle(he).idx() != cut[i - 1].idx());
			}
			he = mesh->next_halfedge_handle(mesh->find_halfedge(cut[cut.size() - 2], cut.back()));
			f = mesh->face_handle(he);
			fid = f.idx();
			while (new_f_flag[fid])
			{
				if (info[fid].empty())
				{
					cutf_flag[fid] = true;
					calc_vector(true);
				}
				he = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(he));
				f = mesh->face_handle(he);
				fid = f.idx();
			}

			for (auto fa : new_face)
			{
				if (cutf_flag[fa.idx()])
					continue;
				he = mesh->fh_begin(fa).handle();
				f = fa;
				fid = f.idx();
				calc_vector(false);
			}
		}
#if PRINT_DEBUG_INFO
		dprint("标记与cut相关的顶点和面，计算装配矩阵需要的数据");
#endif

		int region_vertex_size = region_vertex.size();
		//计算idmap
		{
			int count = region_vertex_size;
			for (auto vv : new_vertex)
				vidmap[vv.idx()] = count++;
		}

		int new_vertex_size = new_vertex.size();
		int new_face_size = new_face.size();
		std::vector<Eigen::Triplet<double>> triple;
		std::vector<double> w(nv, 0);
		std::vector<double> size_ratio(nf, 1.0);
		uv[0].conservativeResize(region_vertex_size + new_vertex_size); uv[0].tail(new_vertex_size).setZero();
		uv[1].conservativeResize(region_vertex_size + new_vertex_size); uv[1].tail(new_vertex_size).setZero();

		//const auto& crossfield = cf->getCrossField();
		for (auto v : new_vertex)
		{
			int vid = v.idx();
			int vm = vidmap[vid];
			for (auto vf = mesh->vf_begin(v); vf != mesh->vf_end(v); ++vf)
			{
				int vf_id = vf->idx();
				if (!new_f_flag[vf_id])
					continue;
				Eigen::Vector3d& R = info[vf_id][v].second;
				for (const auto& f_info : info[vf_id])
				{
					double dot_ = R.dot(f_info.second.second);
					int fvid = f_info.first.idx();
					//w[fvid] += dot_;
					if (f_info.second.first)
						uv[0](vm) -= dot_;
					if (!new_v_flag[fvid])
					{
						uv[0](vm) -= dot_ * GetU(fvid);
						uv[1](vm) -= dot_ * GetV(fvid);
					}
					else
						w[fvid] += dot_;
				}
				//uv[0](vm) += crossfield.col(4 * vf_id).dot(R) * size_ratio[vf_id];
				//uv[1](vm) += crossfield.col(4 * vf_id + 1).dot(R) * size_ratio[vf_id];
				uv[0](vm) += x_axis.col(vf_id).dot(R) * size_ratio[vf_id];
				uv[1](vm) += y_axis.col(vf_id).dot(R) * size_ratio[vf_id];
			}
			triple.emplace_back(vm - region_vertex_size, vm - region_vertex_size, w[vid]);
			w[vid] = 0;
			for (auto vv = mesh->vv_begin(v); vv != mesh->vv_end(v); ++vv)
			{
				int vvid = vv->idx();
				if (!new_v_flag[vvid])
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
#if PRINT_DEBUG_INFO
		dprint("计算参数化");
#endif
	}

#pragma region initialization
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
			if (weight(shift, voh->idx()) < YYSS_INFINITE)
			{
				int toid = mesh->to_vertex_handle(voh.handle()).idx();
				distance[toid] = weight(shift, voh->idx());
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
				double w = weight(shift, voh.idx());// w *= w;
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
			cf = new crossField(field_file);
		}
		else
		{
			file_reader.close();
			cf = new crossField(mesh);
			cf->write_field(field_file);
		}
#else
		cf = new crossField(mesh);
#endif
		dprint("Initialize Field Done!");
	}

	void LoopGen::InitializeGraphWeight(double alpha)
	{
		auto& crossfield = cf->getCrossField();
		auto& matching = cf->getMatching();
		auto& normal = cf->getNormal();
		auto& position = cf->getPosition();

		double doublePI = 2.0 * PI;
		double halfPI = 0.5 * PI;
		double _PI = -1.0 * PI;
		double _halfPI = -0.5 * PI;

		weight.resize(4, mesh->n_halfedges());
		for (auto eitr = mesh->edges_begin(); eitr != mesh->edges_end(); ++eitr)
		{
			auto h0 = mesh->halfedge_handle(eitr.handle(), 0);
			auto h1 = mesh->halfedge_handle(eitr.handle(), 1);
			auto fid = mesh->face_handle(h0).idx();
			auto gid = mesh->face_handle(h1).idx();
			auto& fv = crossfield.col(fid * 4);
			auto& gv = crossfield.col(gid * 4 + matching[h0.idx()]);
			auto ev = position.col(mesh->to_vertex_handle(h0).idx()) - position.col(mesh->from_vertex_handle(h0).idx());
			double arc0 = atan2(ev.cross(fv).dot(normal.col(fid)), ev.dot(fv)); //arc0 += arc0 > 0 ? 0 : doublePI;
			double arc1 = atan2(ev.cross(gv).dot(normal.col(gid)), ev.dot(gv)); //arc1 += arc1 > 0 ? 0 : doublePI;
			double arc = atan2(sin(arc0) + sin(arc1), cos(arc0) + cos(arc1));

			auto& w0 = weight.col(h0.idx());
			auto& w1 = weight.col(h1.idx());
			double s = fabs(sin(arc));
			double c = fabs(cos(arc));

			if (s < c)
			{
				s = sqrt(alpha * s * s + 1);
				c = YYSS_INFINITE;
			}
			else
			{
				s = YYSS_INFINITE;
				c = sqrt(alpha * c * c + 1);
			}
			if (arc >= 0 && arc < halfPI)
				w0 << s, YYSS_INFINITE, YYSS_INFINITE, c;
			else if (arc >= halfPI && arc < PI)
				w0 << YYSS_INFINITE, YYSS_INFINITE, s, c;
			else if (arc >= _PI && arc < _halfPI)
				w0 << YYSS_INFINITE, c, s, YYSS_INFINITE;
			else
				w0 << s, c, YYSS_INFINITE, YYSS_INFINITE;

			switch (matching[h0.idx()])
			{
			case 0:
				w1 << w0(2), w0(3), w0(0), w0(1);
				break;
			case 1:
				w1 << w0(1), w0(2), w0(3), w0(0);
				break;
			case 2:
				w1 << w0(0), w0(1), w0(2), w0(3);
				break;
			case 3:
				w1 << w0(3), w0(0), w0(1), w0(2);
				break;
			}
		}
		dprint("Initialize Graph Weight Done!");
	}

	void LoopGen::InitializePQ()
	{
		timeRecorder tr;
		InfoOnMesh.resize(mesh->n_vertices() * 2);
		for (int i = 0; i < mesh->n_vertices(); ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				InfoOnMesh[2 * i + j].v = mesh->vertex_handle(i);
			}
		}
		dprint(mesh->n_vertices(), InfoOnMesh.size());
		int nedges = mesh->n_edges();

		tr.tog();
		if (COMPUTE_NEW_PLANELOOP || !ReadPlaneLoop(InfoOnMesh, model_name, mesh))
		{
#pragma omp parallel for
			for (int i = 0; i < mesh->n_vertices(); ++i)
			{
				for (int j = 0; j < 2; ++j)
				{
					if (FieldAligned_PlanarLoop(InfoOnMesh[2 * i + j].v, InfoOnMesh[2 * i + j].loop, j))
					{
						RefineLoopByPlanarity(InfoOnMesh[2 * i + j].loop, InfoOnMesh[2 * i + j].pl, j);
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
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(&InfoOnMesh[2 * tovert.idx()], 0));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(&InfoOnMesh[2 * fromvert.idx()], 0));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(&InfoOnMesh[2 * tovert.idx() + 1], 0));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(&InfoOnMesh[2 * fromvert.idx() + 1], 0));
				break;
			case 1:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(&InfoOnMesh[2 * tovert.idx() + 1], 1));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(&InfoOnMesh[2 * fromvert.idx() + 1], 0));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(&InfoOnMesh[2 * tovert.idx()], 0));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(&InfoOnMesh[2 * fromvert.idx()], 1));
				break;
			case 2:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(&InfoOnMesh[2 * tovert.idx()], 1));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(&InfoOnMesh[2 * fromvert.idx()], 1));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(&InfoOnMesh[2 * tovert.idx() + 1], 1));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(&InfoOnMesh[2 * fromvert.idx() + 1], 1));
				break;
			case 3:
				InfoOnMesh[2 * fromvert.idx()].mark.insert(std::make_pair(&InfoOnMesh[2 * tovert.idx() + 1], 0));
				InfoOnMesh[2 * tovert.idx()].mark.insert(std::make_pair(&InfoOnMesh[2 * fromvert.idx() + 1], 1));
				InfoOnMesh[2 * fromvert.idx() + 1].mark.insert(std::make_pair(&InfoOnMesh[2 * tovert.idx()], 1));
				InfoOnMesh[2 * tovert.idx() + 1].mark.insert(std::make_pair(&InfoOnMesh[2 * fromvert.idx()], 0));
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
					auto& tl = fl.mark.find(&InfoOnMesh[2 * tovert.idx() + i]) != fl.mark.end() ?
						InfoOnMesh[2 * tovert.idx() + i] : InfoOnMesh[2 * tovert.idx() + ((i + 1) % 2)];
					if (fl.pl.empty() || tl.pl.empty())
					{
						similarity_energy[2 * k + i] = YYSS_INFINITE;
						continue;
					}

					int flag = fl.mark[&tl];
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
				auto& tl = fl.mark.find(&InfoOnMesh[toid2 + j]) != fl.mark.end() ?
					InfoOnMesh[toid2 + j] : InfoOnMesh[toid2 + ((j + 1) % 2)];
				fl.energy += similarity_energy[i * 2 + j];
				tl.energy += similarity_energy[i * 2 + j];
			}
		}
#pragma omp parallel for
		for (int i = 0; i < InfoOnMesh.size(); ++i)
			InfoOnMesh[i].energy /= mesh->valence(InfoOnMesh[i].v);
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
			if (vv.first->v.idx() == fromid)
			{
				hierarchy_vertex[1].push_back(vv.first);
				break;
			}
		}
		//将plane loop上的点加入hierarchy_vertex中
		for (auto pl_b = pl.begin(); pl_b != pl.end(); ++pl_b)
		{
			auto he = pl_b->h;
			if (mesh->from_vertex_handle(he).idx() == hierarchy_vertex[1].back()->v.idx())
			{
				//hierarchy_vertex[0].push_back(&InfoOnMesh[2*mesh->to_vertex_handle(he).idx()+])
				toid = 2 * mesh->to_vertex_handle(he).idx();
				const auto& mark = hierarchy_vertex[1].back()->mark;
				hierarchy_vertex[0].push_back(mark.find(&InfoOnMesh[toid]) != mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
			}
			else
			{
				fromid = 2 * mesh->from_vertex_handle(he).idx();
				const auto& mark = hierarchy_vertex[0].back()->mark;
				hierarchy_vertex[1].push_back(mark.find(&InfoOnMesh[fromid]) != mark.end() ? &InfoOnMesh[fromid] : &InfoOnMesh[fromid + 1]);
			}
		}
		//将出发点的1邻域加入hierarchy_vertex中
		auto h_b = mesh->next_halfedge_handle(mesh->find_halfedge(hierarchy_vertex[0].back()->v, iov->v));
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
		}
		//advancing_front[0].push_back(std::move(hierarchy_vertex[0]));
		advancing_front[0][0] = std::move(hierarchy_vertex[0]);
		advancing_front[1].push_back(std::move(hierarchy_vertex[1]));

		double energy_threshold = 2.0;

		std::deque<bool> visited_v(mesh->n_vertices(), false);
		visited_v[iov->v.idx()] = true;
		for (auto caf : advancing_front[0].back())
			visited_v[caf->v.idx()] = true;
		for (auto caf : advancing_front[1].back())
			visited_v[caf->v.idx()] = true;

		while (true)
		{
			const auto& current_af = advancing_front[0].back();
			std::vector<InfoOnVertex*> hierarchy;
			for (auto caf : current_af)
			{
				for (const auto& iov_nei : caf->mark)
				{
					int inid = iov_nei.first->v.idx();
					if (visited_v[inid])
						continue;
					if (iov_nei.first->energy > energy_threshold)
						goto target0;
					visited_v[inid] = true;
					hierarchy.push_back(iov_nei.first);
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
					int inid = iov_nei.first->v.idx();
					if (visited_v[inid])
						continue;
					if (iov_nei.first->energy > energy_threshold)
						goto target1;
					visited_v[inid] = true;
					hierarchy.push_back(iov_nei.first);
				}
			}
			advancing_front[1].push_back(std::move(hierarchy));
		}
	target1:;

		auto& new_vertex = lp.GetNewVertex();
		auto& new_face = lp.GetNewFace();
		auto& newv_flag = lp.GetNewVFlag(); newv_flag.resize(nv, false);
		auto& newf_flag = lp.GetNewFFlag(); newf_flag.resize(nf, false);
		auto& regionv_flag = lp.GetRegionVFlag();
		auto& regionf_flag = lp.GetRegionFFlag();
		auto& grow_dir = lp.GetGrowDir();
		for (int i = 0; i < 2; ++i)
		{
			for (auto& ss : advancing_front[i])
			{
				for (auto& tt : ss)
				{
					new_vertex.push_back(tt->v);
					grow_dir[tt->v.idx()] = i;
					newv_flag[tt->v.idx()] = true;
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
		face_tree.push(mesh->face_handle(mesh->voh_begin(iov->v).handle()));
		std::vector<int> ff_id(nf, -1);
		ff_id[face_tree.front().idx()] = iov == &InfoOnMesh[iov->v.idx() * 2] ? 0 : 1;
		auto& matching = cf->getMatching();
		auto& crossfield = cf->getCrossField();
		std::deque<bool> searched(nf, false);
		searched[face_tree.front().idx()] = true;
		
		auto& x_axis = lp.GetXAxis(); 
		auto& y_axis = lp.GetYAxis(); 
		while (!face_tree.empty())
		{
			auto ft = face_tree.front();
			int ftid = ft.idx();
			x_axis.col(ftid) = crossfield.col(4 * ftid + ff_id[ftid]);
			y_axis.col(ftid) = crossfield.col(4 * ftid + (ff_id[ftid] + 1) % 4);
			/*Eigen::Vector3d temp;
			switch (ff_id[ftid])
			{
			case 0:
				break;
			case 1:
				temp = crossfield.col(ftid * 4);
				crossfield.col(ftid * 4) = crossfield.col(ftid * 4 + 1);
				crossfield.col(ftid * 4 + 1) = crossfield.col(ftid * 4 + 2);
				crossfield.col(ftid * 4 + 2) = crossfield.col(ftid * 4 + 3);
				crossfield.col(ftid * 4 + 3) = temp;
				break;
			case 2:
				temp = crossfield.col(ftid * 4);
				crossfield.col(ftid * 4) = crossfield.col(ftid * 4 + 2);
				crossfield.col(ftid * 4 + 2) = temp;
				temp = crossfield.col(ftid * 4 + 1);
				crossfield.col(ftid * 4 + 1) = crossfield.col(ftid * 4 + 3);
				crossfield.col(ftid * 4 + 3) = temp;
				break;
			case 3:
				temp = crossfield.col(ftid * 4 + 3);
				crossfield.col(ftid * 4 + 3) = crossfield.col(ftid * 4 + 2);
				crossfield.col(ftid * 4 + 2) = crossfield.col(ftid * 4 + 1);
				crossfield.col(ftid * 4 + 1) = crossfield.col(ftid * 4);
				crossfield.col(ftid * 4) = temp;
				break;
			}*/
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
		for (auto& pl : lp.GetAllPL()[v.idx()])
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
				double dot_ = fragment.col(i).dot(fragment.col(j));
				if (dot_ > 1.0) dot_ = 1.0;
				if (dot_ < -1.0) dot_ = -1.0;
				sa(count++) = acos(dot_);
			}
		}
		/*for (int i = 0; i < sa.size(); ++i)
			if (fabs(sa(i)) < YYSS_FAIRLY_SMALL)
				sa(i) = YYSS_FAIRLY_SMALL;*/
	}

	bool LoopGen::SpreadSubRegion(LocalParametrization& lp, bool grow_flag[2])
	{
		auto& newv_flag = lp.GetNewVFlag();
		auto& newf_flag = lp.GetNewFFlag();
		auto& new_vertex = lp.GetNewVertex();
		auto& new_face = lp.GetNewFace();
		std::deque<bool> visited_v = lp.GetRegionVFlag();
		std::deque<bool> visited_f = lp.GetRegionFFlag();

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
		if (lp.GetRegionVertex().size() == 1)
			RefineLoopByParametrization(lp.GetRegionVertex().front(), lp, visited_v, visited_f);
		std::vector<VertexHandle> vertex_cache;
		std::deque<bool> vertex_cache_flag = lp.GetRegionVFlag();
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
		if (vertex_cache.empty())
			return false;

		auto& regionv_flag = lp.GetRegionVFlag();
		auto& regionf_flag = lp.GetRegionFFlag();
		auto& grow_dir = lp.GetGrowDir();

#if 0
		double energy_threshold = 0.16; 
		std::ofstream file_writer[2];
		file_writer[0].open("..//resource//statics//" + model_name + ".statics00", std::ios::out | std::ios::app);
		file_writer[1].open("..//resource//statics//" + model_name + ".statics11", std::ios::out | std::ios::app);
		int exceed[2] = { 0,0 };
		std::deque<bool> if_similarity_energy_low(mesh->n_vertices(), true);
		double max_en = -1;
		for (int i = 0; i < vertex_cache.size(); ++i)
		{
			//dprint(i);
			int new_id = vertex_cache[i].idx();
			int grow_id = grow_dir[new_id] ? 1 : 0;
			if (!grow_flag[grow_id])
			{
				if_similarity_energy_low[new_id] = false;
				continue;
			}
			InfoOnVertex* iov0 = &InfoOnMesh[new_id * 2];
			Eigen::VectorXd sa0;
			int loop_fragment_num = iov0->pl.size();
			AssembleSimilarityAngle(iov0, sa0, lp, loop_fragment_num);
			for (auto vv = mesh->vv_begin(vertex_cache[i]); vv != mesh->vv_end(vertex_cache[i]); ++vv)
			{
				int vvid = vv->idx();
				InfoOnVertex* iov1 = &InfoOnMesh[vvid * 2];
				if (!vertex_cache_flag[vvid] || (!regionv_flag[vvid] && new_id > vvid))
					continue;
				Eigen::VectorXd sa1;
				AssembleSimilarityAngle(iov1, sa1, lp, loop_fragment_num);
				double sum = 0;
				for (int j = 0; j < sa0.size(); ++j) {
					//dprint(j, sa0(j), sa1(j), fabs(sin(sa0(j) - sa1(j))));
					//sum += std::min(100.0, fabs(sa0(j) / sa1(j)) + fabs(sa1(j) / sa0(j)) - 2.0);
					sum += fabs(sin(sa0(j) - sa1(j))) + 1 - cos(sa0(j) - sa1(j));
				}
				file_writer[0] << sa0.size() << std::endl;
				file_writer[1] << sum / sa0.size() << std::endl;
				//dprint(sum / sa0.size());
				max_en = std::max(sum / sa0.size(), max_en);
				if (sum > energy_threshold * sa0.size())
				{
					if (exceed[grow_id] > 2)
					{
						if_similarity_energy_low[new_id] = false;
						grow_flag[grow_id] = false;
					}
					++exceed[grow_id];
				}
			}
		}
		file_writer[0].close();
		file_writer[1].close();
		dprint("max ene:", max_en);
#else
		//检测相似性能量
		auto& normal_similarity_angle = lp.GetNormalSimilarityAngle();
		//int loop_fragment_num = InfoOnMesh[2 * lp.GetRegionVertex().front().idx()].pl.size();
		auto& all_pl = lp.GetAllPL();
		int loop_fragment_num = all_pl[lp.GetRegionVertex().front().idx()].size();
		if (vertex_cache_flag[lp.GetRegionVertex().front().idx()] && !lp.has_nsa)
		{
			lp.has_nsa = true;
			AssembleSimilarityAngle(lp.GetRegionVertex().front(), normal_similarity_angle, lp, loop_fragment_num);
		}
		/*Eigen::VectorXd inverse_normal_similarity_angle(normal_similarity_angle.size());
		for (int i = 0; i < normal_similarity_angle.size(); ++i)
			inverse_normal_similarity_angle(i) = 1.0 / normal_similarity_angle(i);*/

		Eigen::VectorXd similarity_angle_sum(normal_similarity_angle.size()); similarity_angle_sum.setZero();
		double energy_threshold = 0.18;
		//static omp_lock_t lock;
		//omp_init_lock(&lock);
//#pragma omp parallel for

		int exceed[2] = { 0,0 };
		std::deque<bool> if_similarity_energy_low(mesh->n_vertices(), true);
		for (int i = 0; i < vertex_cache.size(); ++i)
		{
			//dprint(i);
			int new_id = vertex_cache[i].idx();
			int grow_id = grow_dir[new_id] ? 1 : 0;
			if (!grow_flag[grow_id])
				continue;
			Eigen::VectorXd similarity_angle;
			AssembleSimilarityAngle(vertex_cache[i], similarity_angle, lp, loop_fragment_num);
			double sum = 0;
			for (int j = 0; j < similarity_angle.size(); ++j) {
				//dprint(j);
				//dprint(similarity_angle(j), normal_similarity_angle(j), inverse_normal_similarity_angle(j));
				//sum += std::min(100.0, fabs(normal_similarity_angle(j) / similarity_angle(j)) + fabs(similarity_angle(j) * inverse_normal_similarity_angle(j)) - 2.0);
				sum += fabs(sin(normal_similarity_angle(j) - similarity_angle(j))) + 1 - cos(normal_similarity_angle(j) - similarity_angle(j));
			}
			//dprint(sum / similarity_angle.size());
			if (sum > energy_threshold * similarity_angle.size())
			{
				if (exceed[grow_id] > 2)
				{
					if_similarity_energy_low[new_id] = false;
					grow_flag[grow_id] = false;
				}
				++exceed[grow_id];
			}
			//omp_set_lock(&lock);
			//similarity_angle_sum += similarity_angle;
			//omp_unset_lock(&lock);
		}
		//omp_destroy_lock(&lock);
		//normal_similarity_angle = ( normal_similarity_angle * lp.GetRegionVertex().size() + similarity_angle_sum) / (lp.GetRegionVertex().size() + vertex_cache.size());

#if PRINT_DEBUG_INFO
		dprint("检测相似性能量");
#endif

#endif

		//更新新区域
		int count = lp.GetRegionVertex().size();
		int begin_ = count;
		auto& u_para = lp.GetU();
		auto& v_para = lp.GetV();
		auto& region_vertex = lp.GetRegionVertex();
		auto& vidmap = lp.GetVidMap();
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

		auto& region_face = lp.GetRegionFace();
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
		if (!grow_flag[0] && !grow_flag[1])
			return false;

		//扩展advancing_front
		int extend_layer = 3;
		newv_flag.resize(mesh->n_vertices(), false);
		new_vertex.clear();
		for (int i = begin_; i < count; ++i)
		{
			auto rvi = region_vertex[i];
			bool growid = grow_dir[rvi.idx()];
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
				bool growid = grow_dir[rvj.idx()];
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
		if (!cf->mesh)
			cf->init(mesh);

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
		auto& x_axis = lp.GetXAxis();
		auto& y_axis = lp.GetYAxis();
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
				//e_f *= e_f; e_f *= e_f;
				//e_g *= e_g; e_g *= e_g;
				if (opt_flag[fid])
				{
					triple.emplace_back(count, idmap[fid], e_f);
				}
				else
				{
					//COMPLEX dir = COMPLEX(crossfield.col(fid * 4).dot(faceBase.col(2 * fid)), crossfield.col(fid * 4).dot(faceBase.col(2 * fid + 1)));
					//dir *= dir; dir *= dir;
					COMPLEX dir = COMPLEX(x_axis.col(fid).dot(faceBase.col(2 * fid)), x_axis.col(fid).dot(faceBase.col(2 * fid + 1)));
					b(count) -= e_f * dir;
				}
				if (opt_flag[gid])
				{
					triple.emplace_back(count, idmap[gid], -e_g);
				}
				else
				{
					//COMPLEX dir = COMPLEX(crossfield.col(gid * 4).dot(faceBase.col(2 * gid)), crossfield.col(gid * 4).dot(faceBase.col(2 * gid + 1)));
					//dir *= dir; dir *= dir;
					COMPLEX dir = COMPLEX(x_axis.col(gid).dot(faceBase.col(2 * gid)), x_axis.col(gid).dot(faceBase.col(2 * gid + 1)));
					b(count) = +e_g * dir;
				}
				++count;
			}
		}

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> slu;
		Eigen::SparseMatrix<COMPLEX> A(count, opt_face.size());
		A.setFromTriplets(triple.begin(), triple.end());
		Eigen::SparseMatrix<COMPLEX> AT = A.adjoint();
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> solver;
		solver.compute(AT * A);
		b = solver.solve(AT * b.head(count));
		auto& crossfield = cf->getCrossField();
		for (auto f : opt_face)
		{
			int fid = f.idx();
			double theta = std::arg(b(idmap[fid]));// *0.25;
			//dprint(b(idmap[fid]), theta);
			x_axis.col(fid) = faceBase.col(fid * 2) * cos(theta) + faceBase.col(fid * 2 + 1) * sin(theta);
			y_axis.col(fid) = faceBase.col(fid * 2) * cos(theta + 0.5 * PI) + faceBase.col(fid * 2 + 1) * sin(theta + 0.5 * PI);
		}
	}

	void LoopGen::ConstructRegionCut(InfoOnVertex* iov, std::deque<bool>& visited, std::vector<VertexHandle>& cut)
	{
		cut.clear();
		cut.push_back(iov->v);
		HalfedgeHandle prevhe; prevhe.invalidate();
		const auto& matching = cf->getMatching();
		int shift = iov != &InfoOnMesh[iov->v.idx() * 2] ? 2 : 1;
		for (int s = shift; s < 5; shift += 2, s = shift)
		{
			int smark = s;
			HalfedgeHandle hb = mesh->voh_begin(iov->v);
			int hb_idx = hb.idx();
			while (true)
			{
				double w = YYSS_INFINITE;
				do
				{
					if (weight(s, hb.idx()) < w)
					{
						w = weight(s, hb.idx());
						prevhe = hb;
						smark = s;
					}
					s += matching[hb.idx()]; s %= 4;
					hb = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(hb));
				} while (hb.idx() != hb_idx);
				if (!visited[mesh->to_vertex_handle(prevhe).idx()])
					break;
				cut.push_back(mesh->to_vertex_handle(prevhe));
				//hb = mesh->opposite_halfedge_handle(prevhe);
				hb_idx = mesh->opposite_halfedge_handle(prevhe).idx();
				s = smark;
				hb = mesh->next_halfedge_handle(prevhe);
			}
			std::reverse(cut.begin(), cut.end());
		}
	}

	void LoopGen::OptimizeLoop()
	{
#if 1
		int vid_set[11] = { 16389,31376,14727,1701,6417,15421,38796,16429,20133,22572,5860};
		cf->constraintId.clear();
		cf->constraintVector.resize(0, 0);
		for (int i = 0; i < 11; ++i)
		{
			dprint("iiiiiiiiiiiiiiiiiiiiiiiii", i, vid_set[i]);
			//if (i == 1 || i == 3)
				//continue;
			InfoOnVertex* iov = InfoOnMesh[vid_set[i] * 2].energy < InfoOnMesh[vid_set[i] * 2 + 1].energy ? &InfoOnMesh[vid_set[i] * 2] : &InfoOnMesh[vid_set[i] * 2 + 1];

			//std::vector<std::vector<InfoOnVertex*>> advancing_front[2];

			LocalParametrization lp(*mesh, *cf, iov->v/*, iov == &InfoOnMesh[iov->v.idx() * 2] ? 0 : 1*/);
			ConstructInitialRegion(iov, lp);
			bool grow_flag[2] = { true,true };
			//
			int itertimes = 0;
			int rr = 0;
			do
			{
				std::deque<bool> visited_v = lp.GetRegionVFlag();
				for (auto ver : lp.GetNewVertex())
					visited_v[ver.idx()] = true;
				ConstructRegionCut(iov, visited_v, lp.GetCut());
#if PRINT_DEBUG_INFO
				dprint("计算cut");
#endif
				if (rr == 4 && vid_set[i] == 14727) break;
				++rr;
				dprint(rr);
				lp.run();
			} while (SpreadSubRegion(lp, grow_flag));
			all_plane_loop.insert(all_plane_loop.end(), lp.GetAllPL().begin(), lp.GetAllPL().end());
			region_vertex.insert(region_vertex.end(), lp.GetRegionVertex().begin(), lp.GetRegionVertex().end());
			auto& region_face = lp.GetRegionFace();
			cf->setOuterConstraint(region_face, lp.GetXAxis());
			sub_vertex.insert(sub_vertex.end(), lp.GetNewVertex().begin(), lp.GetNewVertex().end());
			sub_face.insert(sub_face.end(), lp.GetNewFace().begin(), lp.GetNewFace().end());
		}
		cf->runPolynomial();
#else
		int vvvid = 16389;
		InfoOnVertex* iov = InfoOnMesh[vvvid * 2].energy < InfoOnMesh[vvvid * 2 + 1].energy ? &InfoOnMesh[vvvid * 2] : &InfoOnMesh[vvvid * 2 + 1];

		std::vector<std::vector<InfoOnVertex*>> advancing_front[2];

		LocalParametrization lp(*mesh, *cf, iov->v);
		ConstructInitialRegion(iov, lp);
		bool grow_flag[2] = { true,true };
		int shift = iov != &InfoOnMesh[iov->v.idx() * 2] ? 2 : 1;
		int itertimes = 0;
		int rr = 0;
		do
		{
			std::deque<bool> visited_v = lp.GetRegionVFlag();
			for (auto ver : lp.GetNewVertex())
				visited_v[ver.idx()] = true;
			ConstructRegionCut(iov, visited_v, lp.GetCut());
			dprint("计算cut");
			//if (rr == 1) break;
			++rr;
			dprint(rr);
			lp.run();
		} while (SpreadSubRegion(lp, grow_flag));
		all_plane_loop.insert(all_plane_loop.end(), lp.GetAllPL().begin(), lp.GetAllPL().end());
		sub_vertex.insert(sub_vertex.end(), lp.GetNewVertex().begin(), lp.GetNewVertex().end());
		sub_face.insert(sub_face.end(), lp.GetNewFace().begin(), lp.GetNewFace().end());
		region_vertex = lp.GetRegionVertex();
		idmap = lp.GetVidMap();
		uv_para[0] = std::move(lp.uv[0]);
		uv_para[1] = std::move(lp.uv[1]);
#endif
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
		if (id[0] != 1|| id[1] != 1)
		{
			dprint("error in repairing loop:", loop[0].idx());
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
			if (v0(0)*v1[0] + v0(1)*v1[1] + v0(2)*v1[2] < 0)
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
		lp.GetAllPL()[v.idx()] = std::move(planar_loop);
		//iov.pl = std::move(planar_loop);
		return true;
	}

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
		auto loopLength = [&](Eigen::Matrix3Xd& loop, Eigen::VectorXd &seg, int mark)
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
		auto assembleFragment = [&](Eigen::Matrix3Xd& fragment, double u0, int mark, Eigen::Matrix3Xd &loop)
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
				dot = fabs(fragment0.col(i).dot(fragment0.col(j)))/(fabs((fragment1.col(i).dot(fragment1.col(j)))) + YYSS_FAIRLY_SMALL);
				sum += std::min(100.0, dot + 1.0 / (dot + YYSS_FAIRLY_SMALL) - 2);
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
}

////割缝处shift的行
			//double mu = 0;
			//for (auto f : face)
			//{
			//	int fid = f.idx();
			//	if (!cutf_flag[fid])
			//		continue;
			//	Eigen::Vector3d R0(0, 0, 0);
			//	for (const auto& f_info : info[fid])
			//	{
			//		if (f_info.second.first)
			//			R0 += f_info.second.second;
			//	}
			//	mu += R0.squaredNorm();
			//	for (const auto& f_info : info[fid])
			//	{
			//		w[f_info.first.idx()] += R0.dot(f_info.second.second);
			//	}
			//	uv[0](v_size) += R0.dot(crossfield.col(4 * fid + ff_id[fid]));
			//}
			//triple.emplace_back(v_size, v_size, mu);
			//for (int i = 0; i < nv; ++i)
			//{
			//	if (std::fabs(w[i]) < YYSS_FAIRLY_SMALL)
			//		continue;
			//	triple.emplace_back(v_size, vidmap[i], w[i]);
			//}

			////面上比例系数的行
			//for (auto ff : face)
			//{
			//	int ffid = ff.idx();
			//	int ffidmap = fidmap[ffid];
			//	auto& F0 = crossfield.col(4 * ffid + ff_id[ffid]);
			//	triple.emplace_back(ffidmap, ffidmap, F0.dot(F0));
			//	for (const auto& f_info : info[ffid])
			//	{
			//		double dot_ = F0.dot(f_info.second.second);
			//		triple.emplace_back(ffidmap, vidmap[f_info.first.idx()], -dot_);
			//		if (f_info.second.first)
			//		{
			//			uv[0](ffidmap) += dot_;
			//		}
			//	}
			//}

//int v_size = vertex.size();
//int f_size = face.size();
//std::vector<Eigen::Triplet<double>> triple;
//std::vector<double> w(nv, 0);
//std::vector<double> size_ratio(nf, 1.0);
//uv[0].resize(v_size);
//uv[0].setZero();
//const auto& crossfield = cf->getCrossField();
//int vertex_front_id = vertex.front().idx();
//for (int itertimes = 0; itertimes < 1; ++itertimes)
//{
//	//跟顶点相关的行
//	for (auto vv : vertex)
//	{
//		//dprint(vv.idx());
//		int vvid = vv.idx();
//		if (vvid == vertex_front_id)
//			continue;
//		int vvidmap = vidmap[vvid];
//		for (auto vf = mesh->vf_begin(vv); vf != mesh->vf_end(vv); ++vf)
//		{
//			int vf_id = vf->idx();
//			if (!f_flag[vf_id])
//				continue;
//			Eigen::Vector3d& R0 = info[vf->idx()][vv].second;
//			for (const auto& f_info : info[vf_id])
//			{
//				double dot_ = R0.dot(f_info.second.second);
//				w[f_info.first.idx()] += dot_;
//				//triple.emplace_back(vvidmap, fidmap[vf_id], -R0.dot(crossfield.col(4 * vf_id + ff_id[vf_id])));
//				if (f_info.second.first)
//				{
//					uv[0](vvidmap) -= 10*dot_;
//				}
//			}
//			uv[0](vvidmap) += size_ratio[vf_id] * R0.dot(crossfield.col(4 * vf_id + ff_id[vf_id]));
//		}

//		/*if (vvid == vertex.front().idx())
//			w[vvid] += 1.0;*/
//		triple.emplace_back(vvidmap - 1, vvidmap - 1, w[vvid]);
//		w[vvid] = 0;
//		for (auto vvv : mesh->vv_range(vv))
//		{
//			int vvvid = vvv.idx();
//			if (vvvid == vertex_front_id || !v_flag[vvvid])
//				continue;
//			triple.emplace_back(vvidmap - 1, vidmap[vvvid] - 1, w[vvvid]);
//			w[vvvid] = 0;
//		}
//	}
//	
//	Eigen::SparseMatrix<double> A(v_size - 1, v_size - 1);
//	A.setFromTriplets(triple.begin(), triple.end());
//	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//	solver.compute(A);
//	uv[0].tail(v_size - 1) = solver.solve(uv[0].tail(v_size - 1));
//}


		//计算每个面上的场的shift
		//{

		/*std::queue<FaceHandle> face_tree;
		for (auto vf : mesh->vf_range(v))
		{
			if (region_f_flag[vf.idx()])
			{
				face_tree.push(vf);
				break;
			}
		}
		std::deque<bool> search_f(nf, false);
		search_f[face_tree.front().idx()] = true;
		face_field_shift[face_tree.front().idx()] = shift;
		while (!face_tree.empty())
		{
			auto fc = face_tree.front();
			face_tree.pop();
			for (auto fh = mesh->fh_begin(fc); fh != mesh->fh_end(fc); ++fh)
			{
				auto oppo_f = mesh->face_handle(mesh->opposite_halfedge_handle(fh)).idx();
				if (search_f[oppo_f] || !region_f_flag[oppo_f])
					continue;
				search_f[oppo_f] = true;
				face_field_shift[oppo_f] = (face_field_shift[fc.idx()] + matching[fh->idx()]) % 4;
				face_tree.push(mesh->face_handle(oppo_f));
			}
		}*/
		//}


		//IsGood(&InfoOnMesh[9935 * 2], &InfoOnMesh[36412 * 2], lp);
#if 0
		//检测将要加入区域的to_add_vertex的相似性
#pragma omp parallel for
		//for (auto new_v : vertex_cache)
for (int i = 0; i < vertex_cache.size(); ++i)
{
	//dprint(i);
	int new_id = vertex_cache[i].idx();
	int grow_id = grow_dir[new_id] ? 1 : 0;
	if (!grow_flag[grow_id])
		continue;
	InfoOnVertex* iov0 = &InfoOnMesh[new_id * 2];
	for (auto vv = mesh->vv_begin(vertex_cache[i]); vv != mesh->vv_end(vertex_cache[i]); ++vv)
	{
		int vvid = vv->idx();
		InfoOnVertex* iov1 = &InfoOnMesh[vvid * 2];
		if (!vertex_cache_flag[vvid] || (!regionv_flag[vvid] && new_id > vvid))
			continue;
		grow_flag[grow_id] = IsGood(iov0, iov1, lp);
		/*if (!grow_flag[grow_id])
		{
			dprint(new_id, vvid);
			int p = 0;
		}*/
	}
	/*dprint(new_id);
	for (auto vv = mesh->vv_begin(new_v); vv != mesh->vv_end(new_v); ++vv)
	{
		dprint(vv->idx());
		InfoOnVertex* iov1 = &InfoOnMesh[vv->idx() * 2];
		if (new_id > vv->idx() || !visited_v[vv->idx()])
			continue;
		grow_flag[grow_id] = IsGood(iov0, iov1, lp);
		dprint("flag", grow_flag[grow_id]);
	}*/
}
dprint("检测将要加入区域的to_add_vertex的相似性");
dprint(grow_flag[0], grow_flag[1], "grow flag");
//如果没有增长区域，则返回false
if (!grow_flag[0] && !grow_flag[1])
return false;
#endif
