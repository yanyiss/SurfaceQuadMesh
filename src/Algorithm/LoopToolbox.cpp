#include "LoopToolbox.h"
namespace LoopGen
{
	bool LoopGen::FieldAligned_PlanarLoop(VertexHandle v, std::vector<int>& loop, double dist, int shift)
	{
		int nv = m->n_vertices();
		int vid = v.idx();
		bool* visited = new bool[nv];
		double* distance = new double[nv];
		HalfedgeHandle* prev = new HalfedgeHandle[nv];
		for (int i = 0; i < nv; ++i) {
			visited[i] = false;
			distance[i] = DBL_MAX;
		}
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
		auto& normal = cf->getNormal();
		/*
		注意测试
		这里的
		法向
		是否正确！！！！！！！！！！！！！！！！！！！！！！！！！！！！
		*/
		double halfPI = PI * 0.5;
		double doublePI = PI * 2.0;
		double triple_halfPI = halfPI * 3.0;
		for (auto voh = m->voh_begin(v); voh != m->voh_end(v); ++voh)
		{
			int fid = m->face_handle(voh.handle()).idx();
			int gid = m->face_handle(m->opposite_halfedge_handle(voh.handle())).idx();
			plane_normal += crossfield.col(4 * fid + shift + 1);
			auto& v0 = crossfield.col(4 * fid + shift);
			int shift_ = shift;
			shift += matching[voh->idx()];
			shift %= 4;
			auto& v1 = crossfield.col(4 * gid + shift);
			int toid = m->to_vertex_handle(voh.handle()).idx();
			auto ev = (position.col(toid) - position.col(vid)).normalized();

			//double arc = atan2(ev.cross(v0).dot(normal.col(m->face_handle(voh.handle()).idx())) +
			//	ev.cross(v1).dot(normal.col(m->face_handle(m->opposite_halfedge_handle(voh.handle())).idx())), ev.dot(v0 + v1));
			//arc += arc > 0 ? 0 : 2 * PI;
			double arc0 = atan2(ev.cross(v0).dot(normal.col(fid)), ev.dot(v0)); arc0 += arc0 > 0 ? 0 : doublePI;
			double arc1 = atan2(ev.cross(v1).dot(normal.col(gid)), ev.dot(v1)); arc1 += arc1 > 0 ? 0 : doublePI;
			double arc = 0.5 * fabs(arc0 + arc1);
			if (arc >= halfPI && arc <= triple_halfPI)
				continue;

			distance[toid] = fabs(sin(arc));
			pq.emplace(toid, shift_, distance[toid], ++count[toid]);
			prev[toid] = voh.handle();
			visited[toid] = true;
		}
		plane_normal.normalize();

		while (true)
		{
			VertexPQ vert;
			do
			{
				if (pq.empty()) return false;
				vert = pq.top();
				pq.pop();
			} while (vert.count != count[vert.id]);


			int fromid = vert.id;
			if (fromid == vid)
			{
				break;
			}

			shift = vert.shift;
			fromid = vert.id;
			auto voh = m->next_halfedge_handle(prev[fromid]);
			int valence = m->valence(m->vertex_handle(fromid)) - 1;
			for (int i = 0; i < valence; ++i)
			{
				int toid = m->to_vertex_handle(voh).idx();
				if (visited[toid])
					continue;
				int fid = m->face_handle(voh).idx();
				int gid = m->face_handle(m->opposite_halfedge_handle(voh)).idx();
				//plane_normal += crossfield.col(4 * fid + shift + 1);
				auto& v0 = crossfield.col(4 * fid + shift);
				int shift_ = shift;
				shift += matching[voh.idx()];
				shift %= 4;
				auto& v1 = crossfield.col(4 * gid + shift);
				auto ev = (position.col(toid) - position.col(fromid)).normalized();

				//double arc = atan2(ev.cross(v0).dot(normal.col(m->face_handle(voh.handle()).idx())) +
				//	ev.cross(v1).dot(normal.col(m->face_handle(m->opposite_halfedge_handle(voh.handle())).idx())), ev.dot(v0 + v1));
				//arc += arc > 0 ? 0 : 2 * PI;
				double arc0 = atan2(ev.cross(v0).dot(normal.col(fid)), ev.dot(v0)); arc0 += arc0 > 0 ? 0 : doublePI;
				double arc1 = atan2(ev.cross(v1).dot(normal.col(gid)), ev.dot(v1)); arc1 += arc1 > 0 ? 0 : doublePI;
				double arc = 0.5 * fabs(arc0 + arc1);
				if (arc >= halfPI && arc <= triple_halfPI)
					continue;

				double s = sin(arc);
				double c = ev.dot(plane_normal);
				double w = sqrt(s * s + c * c / ev.squaredNorm());
				if (distance[fromid] + w < distance[toid])
				{
					distance[toid] = distance[fromid] + w;
					pq.emplace(toid, shift_, distance[toid], ++count[toid]);
					prev[toid] = voh;
				}
				/*pq.emplace(toid, shift_, distance[toid], ++count[toid]);
				prev[toid] = voh.handle();*/
				visited[fromid] = true;
				voh = m->next_halfedge_handle(m->opposite_halfedge_handle(voh));
			}

			//cfid = faceAlignId[prev[fromid].opp().face().idx()];
			//int valence = mesh->valence(mesh->vertex_handle(fromid)) - 1;
			//auto tvoh = prev[fromid].opp().prev().opp();
			//for (int i = 0; i < valence; ++i)
			//{
			//	cfid += matching[tvoh.idx()]; cfid %= 4;
			//	faceAlignId[tvoh.face().idx()] = cfid;
			//	if (halfedgeAlignId[tvoh.idx()] + cfid == 4)
			//	{
			//		int toid = tvoh.to().idx();
			//		if (!visited[toid] && distance[fromid] + weight(fromid, toid) < distance[toid])
			//		{
			//			distance[toid] = distance[fromid] + weight(fromid, toid);
			//			pq.emplace(toid, distance[toid], ++count[toid]);
			//			prev[toid] = tvoh;
			//			//visited[toid] = true;
			//		}
			//	}
			//}
		}
		return true;
	}

	void LoopGen::InitializeField()
	{
		cf = new crossField(m);
		cf->runPolynomial();
	}

	void LoopGen::InitializePQ()
	{
		//compute_principal_curvature(m, cur[0], cur[1], dir[0], dir[1]);
		for (int i = 0; i < 2; ++i)
		{
			//for (auto &k : cur[i])
				//k = fabs(k);
		}
		for (auto v : m->vertices())
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