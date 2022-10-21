#include "LoopToolbox.h"
namespace LoopGen
{
	bool LoopGen::FieldAligned_PlanarLoop(VertexHandle v, std::vector<int>& loop, int shift)
	{
		int nv = mesh->n_vertices();
		int vid = v.idx();
		std::deque<bool> visited(nv, false);
		std::vector<double> distance(nv, DBL_MAX);
		std::vector<HalfedgeHandle> prev(nv);
		shift %= 2;

		struct VertexPQ
		{
			int id;
			//bool visited;
			int shift;
			HalfedgeHandle prev;
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
			if (weight(shift, voh->idx()) < 2)
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
				double w = weight(shift, voh.idx()); w *= w;
				if (w < 2)
				{
					int toid = mesh->to_vertex_handle(voh).idx();
					double dot_ = (position.col(toid) - position.col(fromid)).normalized().dot(plane_normal); dot_ *= dot_;
					w = sqrt(w + dot_);
					if (distance[fromid] + w < distance[toid])
					{
						distance[toid] = distance[fromid] + w;
						pq.emplace(toid, shift, distance[toid], ++count[toid]);
						prev[toid] = voh;
					}
				}
				voh = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(voh));
				shift += matching[voh.idx()]; shift %= 4;
			}
		}

		loop.push_back(vid);
		int previd = mesh->from_vertex_handle(prev[vid]).idx();
		while (previd != vid)
		{
			loop.push_back(previd);
			previd = mesh->from_vertex_handle(prev[previd]).idx();
		}
		loop.push_back(vid);
		return true;
	}

	void LoopGen::InitializeField()
	{
#if 0
		cf = new crossField(mesh);
#else
		cf = new crossField("..//resource//field//vase.field");
#endif
	}

	void LoopGen::InitializeGraphWeight()
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
			if (arc >= 0 && arc < halfPI)
			{
				w0 << s, DBL_MAX, DBL_MAX, c;
				w1 << DBL_MAX, c, s, DBL_MAX;
			}
			else if (arc >= halfPI && arc < PI)
			{
				w0 << DBL_MAX, DBL_MAX, s, c;
				w1 << s, c, DBL_MAX, DBL_MAX;
			}
			else if (arc >= _PI && arc < _halfPI)
			{
				w0 << DBL_MAX, c, s, DBL_MAX;
				w1 << s, DBL_MAX, DBL_MAX, c;
			}
			else
			{
				w0 << s, c, DBL_MAX, DBL_MAX;
				w1 << DBL_MAX, DBL_MAX, s, c;
			}
		}
	}

	void LoopGen::InitializePQ()
	{
		//compute_principal_curvature(m, cur[0], cur[1], dir[0], dir[1]);
		for (int i = 0; i < 2; ++i)
		{
			//for (auto &k : cur[i])
				//k = fabs(k);
		}
		for (auto v : mesh->vertices())
		{
			auto vid = v.idx();
			EnergyOnVertex ev;
			ev.v = v;
			//if (cur[0][vid] < 1.0e-6 || cur[1][vid] < 1.0e-6)
			{
				ev.energy = DBL_MAX;
				pq.push(ev);
				continue;
			}


			//FieldAligned_PlanarLoop(v, dir[cur[0][vid] >= cur[1][vid] ? 0 : 1][vid]);
			LocalParametrization lp;
			ev.energy = 0;
			pq.push(ev);
		}
	}
}