#include "LoopToolbox.h"
#include <Eigen\Cholesky>
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
		std::vector<int> ff_id(mesh->n_faces());
		{
			std::queue<FaceHandle> face_tree;
			face_tree.push(mesh->voh_begin(v).handle().face());
			std::deque<bool> visited_f(mesh->n_faces(), false);
			visited_f[face_tree.front().idx()] = true;
			ff_id[face_tree.front().idx()] = shift;
			while (!face_tree.empty())
			{
				auto f = face_tree.front();
				face_tree.pop();
				for (auto fh = mesh->fh_begin(f); fh != mesh->fh_end(f); ++fh)
				{
					auto oppo_f = mesh->face_handle(mesh->opposite_halfedge_handle(fh)).idx();
					if (visited_f[oppo_f])
						continue;
					visited_f[oppo_f] = true;
					ff_id[oppo_f] = (ff_id[f.idx()] + matching[fh->idx()]) % 4;
					face_tree.push(mesh->face_handle(oppo_f));
				}
			}
		}

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

			//loop.push_back(vert.idmap); loop.push_back(mesh->from_vertex_handle(prev[vert.idmap]).idx());
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
						eov[i] = std::min(eov[i], RefineLoop(InfoOnMesh[2 * i + j].loop, InfoOnMesh[2 * i + j].pl, j));
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

	void LoopGen::ConstructSubRegion(InfoOnVertex* iov, std::vector<std::vector<InfoOnVertex*>> advancing_front[2], std::deque<bool>& visited_v)
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
		visited_v.resize(mesh->n_vertices(), false);

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
		}
	}

	void LoopGen::OptimizeLoop()
	{
		InfoOnVertex* iov = InfoOnMesh[33233 * 2].energy < InfoOnMesh[33233 * 2 + 1].energy ? &InfoOnMesh[33233 * 2] : &InfoOnMesh[33233 * 2 + 1];
		std::vector<std::vector<InfoOnVertex*>> advancing_front[2];
		std::deque<bool> visited_v;
		ConstructSubRegion(iov, advancing_front, visited_v);

		std::deque<bool> visited_f(mesh->n_faces(), false);
		LocalParametrization lp(*mesh, *cf);
		auto& subvertex = lp.GetVertex(); subvertex.clear();
		auto& subface = lp.GetFace(); subface.clear();
		for (auto& ss : advancing_front)
		{
			for (auto& tt : ss)
			{
				for (auto& rr : tt)
				{
					subvertex.push_back(rr->v);
					for (auto vf = mesh->vf_begin(rr->v); vf != mesh->vf_end(rr->v); ++vf)
					{
						if (visited_f[vf->idx()])
							continue;
						visited_f[vf->idx()] = true;
						subface.push_back(vf.handle());
					}
				}
			}
		}
		ConstructRegionCut(iov->v, iov != &InfoOnMesh[iov->v.idx() * 2] ? 0 : 1, visited_v, lp.GetCut());

		sub_vertex = lp.GetVertex();
		sub_face = lp.GetFace();
		sub_cut = lp.GetCut();

		lp.run(iov->v, iov == &InfoOnMesh[iov->v.idx() * 2] ? 0 : 1);
	}

	double LoopGen::RefineLoop(std::vector<VertexHandle>& loop, PlaneLoop& planar_loop, int shift)
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
				/*h[idmap] = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(hitr.handle()));
				++idmap;*/
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