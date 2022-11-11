#include "LoopToolbox.h"
#include <omp.h>
#define YYSS_INFINITE 1.0e12
#define YYSS_FAIRLY_SMALL 1.0e-6
#define COMPUTE_NEW_PLANELOOP 0
#define COMPUTE_NEW_ENERGY 0
namespace LoopGen
{
	void LocalParametrization::run(VertexHandle v, int shift)
	{
		const auto& matching = cf->getMatching();
		int nf = mesh->n_faces();
		//计算每个面上的场的shift
		std::vector<int> ff_id(nf, -1);
		{
			std::queue<FaceHandle> face_tree;
			for (auto vf : mesh->vf_range(v))
			{
				if (f_flag[vf.idx()])
				{
					face_tree.push(vf);
					break;
				}
			}
			std::deque<bool> search_f(nf, false);
			search_f[face_tree.front().idx()] = true;
			ff_id[face_tree.front().idx()] = shift;
			while (!face_tree.empty())
			{
				auto fc = face_tree.front();
				face_tree.pop();
				for (auto fh = mesh->fh_begin(fc); fh != mesh->fh_end(fc); ++fh)
				{
					auto oppo_f = mesh->face_handle(mesh->opposite_halfedge_handle(fh)).idx();
					if (search_f[oppo_f] || !f_flag[oppo_f])
						continue;
					search_f[oppo_f] = true;
					ff_id[oppo_f] = (ff_id[fc.idx()] + matching[fh->idx()]) % 4;
					face_tree.push(mesh->face_handle(oppo_f));
				}
			}
		}

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
				info[fid].insert(std::make_pair(vi, std::make_pair(flag && cutv_flag[vi.idx()], inv_area * normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])))));
				vi = mesh->to_vertex_handle(mesh->next_halfedge_handle(he));
				ev = mesh->calc_edge_vector(he);
				info[fid].insert(std::make_pair(vi, std::make_pair(flag && cutv_flag[vi.idx()], inv_area * normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])))));
				vi = mesh->from_vertex_handle(he);
				ev = mesh->calc_edge_vector(mesh->next_halfedge_handle(he));
				info[fid].insert(std::make_pair(vi, std::make_pair(flag && cutv_flag[vi.idx()], inv_area * normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])))));
			};

			while (f_flag[fid])
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
					if (info[fid].empty())
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
			while (f_flag[fid])
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

			for (auto fa : face)
			{
				if (cutf_flag[fa.idx()])
					continue;
				he = mesh->fh_begin(fa).handle();
				f = fa;
				fid = f.idx();
				calc_vector(false);
			}
		}


		//计算idmap
		vidmap.resize(nv);
		{
			int count = 0;
			for (auto vv : vertex)
				vidmap[vv.idx()] = count++;
			/*for (auto c : cut)
				vidmap[c.idx()] = count++;
			for (auto vv : vertex)
			{
				if (cutv_flag[vv.idx()])
					continue;
				vidmap[vv.idx()] = count++;
			}*/
		}

		int v_size = vertex.size();
		int f_size = face.size();
		std::vector<Eigen::Triplet<double>> triple;
		std::vector<double> w(nv, 0);
		std::vector<double> size_ratio(nf, 1.0);
		uv[0].resize(v_size); uv[0].setZero();
		uv[1].resize(v_size); uv[1].setZero();
		const auto& crossfield = cf->getCrossField();
		int vertex_front_id = vertex.front().idx();

		for (auto vv : vertex)
		{
			//dprint(vv.idx());
			int vvid = vv.idx();
			if (vvid == vertex_front_id)
				continue;
			int vvidmap = vidmap[vvid];
			for (auto vf = mesh->vf_begin(vv); vf != mesh->vf_end(vv); ++vf)
			{
				int vf_id = vf->idx();
				if (!f_flag[vf_id])
					continue;
				Eigen::Vector3d& R0 = info[vf->idx()][vv].second;
				for (const auto& f_info : info[vf_id])
				{
					double dot_ = R0.dot(f_info.second.second);
					w[f_info.first.idx()] += dot_;
					//triple.emplace_back(vvidmap, fidmap[vf_id], -R0.dot(crossfield.col(4 * vf_id + ff_id[vf_id])));
					if (f_info.second.first)
					{
						uv[0](vvidmap) -= dot_;
					}
				}
				uv[0](vvidmap) += size_ratio[vf_id] * R0.dot(crossfield.col(4 * vf_id + ff_id[vf_id]));
				uv[1](vvidmap) += size_ratio[vf_id] * R0.dot(crossfield.col(4 * vf_id + (ff_id[vf_id] + 1) % 4));
			}

			/*if (vvid == vertex.front().idx())
				w[vvid] += 1.0;*/
			triple.emplace_back(vvidmap - 1, vvidmap - 1, w[vvid]);
			w[vvid] = 0;
			for (auto vvv : mesh->vv_range(vv))
			{
				int vvvid = vvv.idx();
				if (vvvid == vertex_front_id || !v_flag[vvvid])
					continue;
				triple.emplace_back(vvidmap - 1, vidmap[vvvid] - 1, w[vvvid]);
				w[vvvid] = 0;
			}
		}
		Eigen::SparseMatrix<double> A(v_size - 1, v_size - 1);
		A.setFromTriplets(triple.begin(), triple.end());
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		uv[0].tail(v_size - 1) = solver.solve(uv[0].tail(v_size - 1));
		uv[1].tail(v_size - 1) = solver.solve(uv[1].tail(v_size - 1));
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

	/*void LoopGen::InitializeAABBTREE()
	{
		ClosestPointSearch::Triangles primitives;
		primitives.reserve(mesh->n_faces());
		OpenMesh::Vec3d p0, p1, p2;
		for (auto fitr = mesh->faces_begin(); fitr != mesh->faces_end(); ++fitr)
		{
			auto fv_iter = mesh->cfv_begin(fitr.handle());
			p0 = mesh->point(fv_iter.handle()); ++fv_iter;
			p1 = mesh->point(fv_iter.handle()); ++fv_iter;
			p2 = mesh->point(fv_iter.handle());
			primitives.emplace_back(Vec3d(p0.data()), Vec3d(p1.data()), Vec3d(p2.data()), fitr.handle());
		}
		aabbtree = new ClosestPointSearch::AABBTree(primitives.begin(), primitives.end());
	}*/

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

			/*if (eitr->idx() == 6097)
			{
				dprint("eitr0");
				dprint(fv.dot(ev.normalized()), gv.dot(ev.normalized()), crossfield.col(4 * gid).dot(ev.normalized()));
				dprint(mesh->from_vertex_handle(h0).idx(), mesh->to_vertex_handle(h0).idx());
				dprint(matching[h0.idx()]);
				int p = 0;
			}*/

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
		eov.resize(mesh->n_vertices(), YYSS_INFINITE);
		timeRecorder tr;
		InfoOnMesh.resize(mesh->n_vertices() * 2);
		int nedges = mesh->n_edges();

		tr.tog();
		if (COMPUTE_NEW_PLANELOOP || !ReadPlaneLoop(InfoOnMesh, model_name, mesh))
		{
#pragma omp parallel for
			for (int i = 0; i < mesh->n_vertices(); ++i)
			{
				for (int j = 0; j < 2; ++j)
				{
					InfoOnMesh[2 * i + j].v = mesh->vertex_handle(i);
					if (FieldAligned_PlanarLoop(InfoOnMesh[2 * i + j].v, InfoOnMesh[2 * i + j].loop, j))
					{
						eov[i] = std::min(eov[i], RefineLoopByPlanarity(InfoOnMesh[2 * i + j].loop, InfoOnMesh[2 * i + j].pl, j));
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
					double e = EvaluateSimilarity(loop0, loop1, u0, id);
					similarity_energy[2 * k + i] = e;
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
		tr.out("time of setting vertex energy");
	}

	void LoopGen::ConstructSubRegion(InfoOnVertex* iov, std::vector<std::vector<InfoOnVertex*>> advancing_front[2])
	{
		auto& pl = iov->pl;
		std::vector<InfoOnVertex*> IOV; 
		IOV.push_back(iov);
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
		hierarchy_vertex[0].shrink_to_fit(); advancing_front[0].push_back(std::move(hierarchy_vertex[0]));
		hierarchy_vertex[1].shrink_to_fit(); advancing_front[1].push_back(std::move(hierarchy_vertex[1]));

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
			hierarchy.reserve(current_af.capacity());
			for (auto caf : current_af)
			{
				const auto& mark = caf->mark;
				for (auto vitr = mesh->vv_begin(caf->v); vitr != mesh->vv_end(caf->v); ++vitr)
				{
					toid = vitr->idx();
					if (visited_v[toid])
						continue;
					visited_v[toid] = true;
					toid *= 2;
					hierarchy.push_back(mark.find(&InfoOnMesh[toid]) != mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
					if (hierarchy.back()->energy > energy_threshold)
						goto target0;
				}
			}
			advancing_front[0].push_back(std::move(hierarchy));
		}
	target0:;
		while (true)
		{
			const auto& current_af = advancing_front[1].back();
			std::vector<InfoOnVertex*> hierarchy;
			hierarchy.reserve(current_af.capacity());
			for (auto caf : current_af)
			{
				const auto& mark = caf->mark;
				for (auto vitr = mesh->vv_begin(caf->v); vitr != mesh->vv_end(caf->v); ++vitr)
				{
					toid = vitr->idx();
					if (visited_v[toid])
						continue;
					visited_v[toid] = true;
					toid *= 2;
					hierarchy.push_back(mark.find(&InfoOnMesh[toid]) != mark.end() ? &InfoOnMesh[toid] : &InfoOnMesh[toid + 1]);
					if (hierarchy.back()->energy > energy_threshold)
						goto target1;
				}
			}
			advancing_front[1].push_back(std::move(hierarchy));
		}
	target1:;
	}

	void LoopGen::ConstructRegionCut(VertexHandle v, int shift, std::deque<bool>& visited, std::vector<VertexHandle>& cut)
	{
		cut.clear();
		cut.push_back(v);
		HalfedgeHandle prevhe; prevhe.invalidate();
		const auto& matching = cf->getMatching();

		for (int s = shift; s < 4; shift += 2, s = shift)
		{
			int smark = s;
			HalfedgeHandle hb = mesh->voh_begin(v);
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
		InfoOnVertex* iov = InfoOnMesh[33233 * 2].energy < InfoOnMesh[33233 * 2 + 1].energy ? &InfoOnMesh[33233 * 2] : &InfoOnMesh[33233 * 2 + 1];
		ConstructSubRegion(iov, advancing_front);

		LocalParametrization lp(*mesh, *cf);
		auto& subvertex = lp.GetVertex(); subvertex.clear();
		auto& subface = lp.GetFace(); subface.clear();
		auto& visited_f = lp.GetFFlag(); visited_f.resize(mesh->n_faces(), false);
		auto& visited_v = lp.GetVFalg(); visited_v.resize(mesh->n_vertices(), false);
		for (auto& ss : advancing_front)
		{
			for (auto& tt : ss)
			{
				for (auto& rr : tt)
				{
					subvertex.push_back(rr->v);
					visited_v[rr->v.idx()] = true;
				}
			}
		}
		for (auto& ss : advancing_front)
		{
			for (auto& tt : ss)
			{
				for (auto& rr : tt)
				{
					for (auto vf = mesh->vf_begin(rr->v); vf != mesh->vf_end(rr->v); ++vf)
					{
						if (visited_f[vf->idx()])
							continue;
						for (auto vfv = mesh->fv_begin(vf.handle()); vfv != mesh->fv_end(vf.handle()); ++vfv)
						{
							if (!visited_v[vfv->idx()])
								goto target;
						}
						visited_f[vf->idx()] = true;
						subface.push_back(vf.handle());
					target:;
					}
				}
			}
		}
		ConstructRegionCut(iov->v, iov != &InfoOnMesh[iov->v.idx() * 2] ? 0 : 1, visited_v, lp.GetCut());

		sub_vertex = lp.GetVertex();
		sub_face = lp.GetFace();
		sub_cut = lp.GetCut();
		
		lp.run(iov->v, iov == &InfoOnMesh[iov->v.idx() * 2] ? 0 : 1);


		uv_para[0] = lp.uv[0];
		uv_para[1] = lp.uv[1];

		std::deque<bool> update_mark(mesh->n_vertices(), false);
		std::queue<InfoOnVertex*> iov_tree;
		iov_tree.push(iov);
		update_mark[iov->v.idx()] = true;
		while (!iov_tree.empty())
		{
			InfoOnVertex* iov_front = iov_tree.front();
			iov_tree.pop();
			RefineLoopByParametrization(*iov_front, lp);
			for (const auto& one_ring : iov_front->mark)
			{
				int orfvid = one_ring.first->v.idx();
				if (update_mark[orfvid] || !visited_v[orfvid])
					continue;
				update_mark[orfvid] = true;
				iov_tree.push(one_ring.first);
			}
		}

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

	void LoopGen::RefineLoopByParametrization(InfoOnVertex& iov, LocalParametrization& lp)
	{
		iov.pl.clear();
		double v_para = lp.GetV(iov.v.idx());
		auto hitr = mesh->voh_begin(iov.v);
		double s0 = v_para - lp.GetV(mesh->to_vertex_handle(mesh->next_halfedge_handle(hitr.handle())).idx());
		double distance[2];
		PointOnHalfedge poh[2];
		int id[2] = { 0,0 };
		for (; hitr != mesh->voh_end(iov.v); ++hitr)
		{
			if (!lp.GetFFlag()[hitr.handle().face().idx()])
				continue;
			double s1 = v_para - lp.GetV(mesh->to_vertex_handle(hitr.handle()).idx());
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
		}
		if (id[0] != 1 || id[1] != 1)
			return;

		////检查出发点到poh[0]的方向与u参数增加的方向是否相同，若不相同，则调换poh[0]和poh[1]
		//auto hb = poh[0].h;
		//const auto& nor = cf->getNormal().col(mesh->face_handle(hb).idx());
		//OpenMesh::Vec3d normal_ = Vec3d(nor(0), nor(1), nor(2));
		//OpenMesh::Vec3d v0(0, 0, 0);
		//for (int i = 0; i < 3; ++i)
		//{
		//	int vid = mesh->to_vertex_handle(mesh->next_halfedge_handle(hb)).idx();
		//	v0 += (lp.GetU(vid) + (lp.GetCutV_Flag()[vid] && lp.GetCutF_Flag()[mesh->face_handle(hb).idx()] ? 1.0 : 0.0))
		//		* normal_.cross(mesh->calc_edge_vector(hb));
		//}
		//auto v1 = poh[0].c * mesh->point(mesh->from_vertex_handle(poh[0].h)) +
		//	(1 - poh[0].c) * mesh->point(mesh->to_vertex_handle(poh[0].h)) - mesh->point(iov.v);
		//if (v0.dot(v1) < 0)
		//{
		//	std::swap(poh[0], poh[1]);
		//	for (int i = 0; i < 2; ++i)
		//	{
		//		poh[i].h = mesh->opposite_halfedge_handle(poh[i].h);
		//		poh[i].c = 1 - poh[i].c;
		//	}
		//	distance[0] = v_para - lp.GetV(mesh->to_vertex_handle(poh[0].h).idx());
		//	distance[1] = v_para - lp.GetV(mesh->from_vertex_handle(poh[0].h).idx());
		//}

		PlaneLoop planar_loop;
		planar_loop.push_back(poh[0]);
		auto h = poh[0].h;
		while (h.idx() != poh[1].h.idx())
		{
			int vid = mesh->from_vertex_handle(mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(h))).idx();
			if (!lp.v_flag[vid])
				//return;
				break;
			double s = v_para - lp.GetV(vid);
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
		iov.pl = std::move(planar_loop);
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

	double LoopGen::ComputeAdjVertexSimilarity(InfoOnVertex& iov0, InfoOnVertex& iov1)
	{
		//首先确定iov0.loop[0]和iov1.loop[0], iov1.loop[1]中哪一个对应, 以及loop的行进方向
		auto& matching = cf->getMatching();
		int index = 0;
		auto h = mesh->find_halfedge(iov0.v, iov1.v);
		auto ht = mesh->voh_begin(iov0.v).handle();
		while (ht.idx() != h.idx())
		{
			ht = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(ht));
			index += matching[ht.idx()];
		}
		h = mesh->next_halfedge_handle(h);
		ht = mesh->voh_begin(iov1.v).handle();
		while (ht.idx() != h.idx())
		{
			ht = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(ht));
			index += 4 - matching[ht.idx()];
		}
		index %= 4;
		//index=
		//0  iov0

		double e0, e1;
		for (int i = 0; i < 2; ++i)
		{

		}
		return std::min(e0, e1);
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

	void LeastSquarePlane(Eigen::VectorXd xyz[3], double plane[4])
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

	double EvaluatePlanarity(Eigen::VectorXd xyz[3], double plane[4])
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

	void boundingXY(Eigen::Matrix3Xd& position, double bounding[4])
	{
		bounding[0] = YYSS_INFINITE; bounding[1] = -YYSS_INFINITE;
		bounding[2] = YYSS_INFINITE; bounding[3] = -YYSS_INFINITE;
		for (int i = 0; i < position.cols(); ++i)
		{
			double r = position(0, i);
			bounding[0] = std::min(bounding[0], r);
			bounding[1] = std::max(bounding[1], r);
			r = position(1, i);
			bounding[2] = std::min(bounding[2], r);
			bounding[3] = std::max(bounding[3], r);
		}
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