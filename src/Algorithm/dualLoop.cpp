#include "dualLoop.h"
namespace QuadLayout
{
#if 0
	void dualLoop::computeGraphWeight(double alpha = 900)
	{
		int nv = mesh->n_vertices();
		double invHalfPi = 2.0 / PI;
		auto &crossfield = cf->getCrossField();
		auto &matching = cf->getMatching();
		auto &position = cf->getPosition();
		Eigen::Matrix3Xd normal(3, mesh->n_faces());
		for (auto &tf : mesh->faces())
		{
			std::vector<int> id; id.reserve(3);
			for (auto &tfv : mesh->fv_range(tf)) id.push_back(tfv.idx());
			normal.col(tf.idx()) = (position.col(id[1]) - position.col(id[0])).cross(position.col(id[2]) - position.col(id[1])).normalized();
		}

		/*for (int i = 0; i < 2; ++i)
		{
			weight[i].resize(nv, nv);
			weight[i].setConstant(DBL_MAX);
		}*/
		weight.resize(nv, nv);
		weight.setConstant(DBL_MAX);
		halfedgeAlignId.resize(mesh->n_halfedges());

		//对于一条半边，只有一个方向使之在图中权重有限，标记为alignId中元素，其他默认为无穷大
		for (auto &te : mesh->edges())
		{
			if (mesh->is_boundary(te))
				continue;
			auto fid = te.h0().face().idx();
			auto gid = te.h1().face().idx();

			auto &fv = crossfield.col(fid * 4);
			auto &gv = crossfield.col(gid * 4 + matching[te.h1().idx()]);
			auto ev = position.col(te.h0().to().idx()) - position.col(te.h0().from().idx());
			double arc = atan2(ev.cross(fv).dot(normal.col(fid)) + ev.cross(gv).dot(normal.col(gid)), ev.dot(fv + gv));
			arc += arc > 0 ? 0 : 2 * PI;
			
			int ind = static_cast<int>((arc + PI * 0.25) * invHalfPi) % 4;
			halfedgeAlignId[te.h0().idx()] = ind;
			halfedgeAlignId[te.h1().idx()] = (ind + 2) % 4;
			arc -= ind * PI * 0.5;

			/*double c = cos(arc); c *= c;
			double s = sin(arc); s *= s;
			double w[2] = { sqrt(c + alpha * s), sqrt(s + alpha * c) };
			for (int i = 0; i < 2; ++i)
			{
				weight[i](te.h0().from().idx(), te.h0().to().idx()) = w[i];
				weight[i](te.h0().to().idx(), te.h0().from().idx()) = w[i];
			}*/
			double w = sqrt(cos(arc)*cos(arc) + alpha * sin(arc)*sin(arc));
			weight(te.h0().from().idx(), te.h0().to().idx()) = w;
			weight(te.h0().to().idx(), te.h0().from().idx()) = w;
		}

	}

	//假设网格没有边界
	bool dualLoop::DijkstraLoop(int vid, int cfid, std::vector<int> &loop, double length)
	{
		recoverDijkstraMark();
		int nv = mesh->n_vertices();
		std::vector<int> count(nv, 0);
		std::vector<int> faceAlignId(mesh->n_faces(), 0);
		auto &matching = cf->getMatching();
		std::priority_queue<VertexPQ, vector<VertexPQ>, std::greater<VertexPQ>> pq;

		//halfedgeAlignId: 从当前半边方向转到所在面第0个crossfield
		//faceAlignId: 从当前面第0个crossfield转到给定场的方向
		for (auto &tvoh : mesh->voh_range(mesh->vertex_handle(vid)))
		{
			faceAlignId[tvoh.face().idx()] = cfid;
			if (cfid + halfedgeAlignId[tvoh.idx()] == 4)
			{
				int toid = tvoh.to().idx();
				pq.emplace(toid, weight(vid, toid), ++count[toid]);
				dprint(count[toid], pq.top().count);
				distance[toid] = weight(vid, toid);
				prev[toid] = tvoh;
			}
			cfid += matching[tvoh.prev().idx()]; cfid %= 4;
		}
		visited[vid] = true;
		while (true)
		{
			VertexPQ tv;
			do
			{
				if (pq.empty()) return false;
				tv = pq.top();
				dprint(tv.count);
				pq.pop();
			} while (tv.count != count[tv.id]);
			int fromid = tv.id;
			if (fromid == vid)
			{
				break;
			}
			cfid = faceAlignId[prev[fromid].opp().face().idx()];
			int valence = mesh->valence(mesh->vertex_handle(fromid)) - 1;
			auto tvoh = prev[fromid].opp().prev().opp();
			for (int i = 0; i < valence; ++i)
			{
				cfid += matching[tvoh.idx()]; cfid %= 4;
				faceAlignId[tvoh.face().idx()] = cfid;
				if (halfedgeAlignId[tvoh.idx()] + cfid == 4)
				{
					int toid = tvoh.to().idx();
					if (!visited[toid] && distance[fromid] + weight(fromid, toid) < distance[toid])
					{
						distance[toid] = distance[fromid] + weight(fromid, toid);
						pq.emplace(toid, distance[toid], ++count[toid]);
						prev[toid] = tvoh;
						//visited[toid] = true;
					}
				}
			}
		}
		
		loop.push_back(vid);
		int previd = prev[vid].from().idx();
		while(previd != vid)
		{
			loop.push_back(previd);
			previd = prev[previd].from().idx();
		}
		length = distance[vid];
		return true;
	}
#endif
	//bool dualLoop::DijkstraLoop(int vid, int cfid, std::vector<int> &path, double length)
	//{
	//	int nv = mesh->n_vertices();
	//	std::vector<int> count(nv, 0);
	//	std::vector<int> hcfid(mesh->n_halfedges(), 0);
	//	auto &matching = cf->getMatching();
	//
	//	std::priority_queue<VertexPQ, vector<VertexPQ>, std::greater<VertexPQ>> pq;
	//	int turnId = cfid;
	//	//在1邻域中挑选与场所在方向在45度以内的边，若没有，才加入其他边
	//	bool hasAlignedEdge = false;
	//	for (auto &tvoh : mesh->voh_range(mesh->vertex_handle(vid)))
	//	{
	//		if (turnId == alignId[tvoh.idx()])
	//		{
	//			int toid = tvoh.to().idx();
	//			++count[toid];
	//			pq.emplace(toid, weight[0](tvoh.idx()), count[toid]);
	//			prev[toid] = tvoh;
	//			distance[toid] = weight[0](tvoh.idx());
	//			hasAlignedEdge = true;
	//			hcfid[tvoh.idx()] = turnId;
	//		}
	//		turnId -= matching[tvoh.prev().idx()]; turnId %= 4;
	//	}
	//	turnId = cfid;
	//	if (!hasAlignedEdge)
	//	{
	//		for (auto &tvoh : mesh->voh_range(mesh->vertex_handle(vid)))
	//		{
	//			if (std::abs(turnId - alignId[tvoh.idx()]) <= 1)
	//			{
	//				int toid = tvoh.to().idx();
	//				++count[toid];
	//				pq.emplace(toid, weight[1](tvoh.idx()), count[toid]);
	//				prev[toid] = tvoh;
	//				distance[toid] = weight[1](tvoh.idx());
	//			}
	//			turnId -= matching[tvoh.prev().idx()]; turnId %= 4;
	//		}
	//	}
	//	visited[vid] = true;
	//	distance[vid] = 0;
	//	while (true)
	//	{
	//		VertexPQ tv;
	//		do
	//		{
	//			if (pq.empty())
	//			{
	//				dprint("fail to find the loop");
	//				return false;
	//			}
	//			tv = pq.top();
	//			pq.pop();
	//		} while (tv.count != count[tv.id]);
	//		int m = tv.id;
	//		visited[m] = true;
	//		hasAlignedEdge = false;
	//		cfid = (hcfid[prev[tv.id].idx()] + matching[prev[tv.id].idx()]) % 4;
	//		turnId = cfid;
	//		int valence = mesh->valence(mesh->vertex_handle(m)) - 1;
	//		for (int i = 0; i < valence; ++i)
	//		{
	//		}
	//		for (auto &tvoh : mesh->voh_range(mesh->vertex_handle(m)))
	//		{
	//			if (tvoh == prev[tv.id])
	//				continue;
	//			int toid = tvoh.to().idx();
	//			if (turnId == alignId[tvoh.idx()] && !visited[toid]  )
	//			{
	//				++count[toid];
	//				pq.emplace(toid, weight[0](tvoh.idx()), count[toid]);
	//				prev[toid] = tvoh;
	//				distance[toid] = weight[0](tvoh.idx());
	//				hasAlignedEdge = true;
	//				hcfid[tvoh.idx()] = turnId;
	//			}
	//			turnId -= matching[tvoh.prev().idx()]; turnId %= 4;
	//		}
	//		turnId = cfid;
	//		if (!hasAlignedEdge)
	//		{
	//			for (auto &tvoh : mesh->voh_range(mesh->vertex_handle(m)))
	//			{
	//				if (std::abs(turnId - alignId[tvoh.idx()]) <= 1)
	//				{
	//					int toid = tvoh.to().idx();
	//					++count[toid];
	//					pq.emplace(toid, weight[1](tvoh.idx()), count[toid]);
	//					prev[toid] = tvoh;
	//					distance[toid] = weight[1](tvoh.idx());
	//					hasAlignedEdge = true;
	//				}
	//				turnId -= matching[tvoh.prev().idx()]; turnId %= 4;
	//			}
	//		}
	//		/*auto it = mesh.vertex_handle(min_id);
	//		auto m = min_id;
	//		for (auto ip = mesh.cvv_begin(it); ip != mesh.cvv_end(it); ++ip)
	//		{
	//			auto n = ip->idx();
	//			if (!visit[n] && dist[m] + admatrix[m][n] < dist[n])
	//			{
	//				dist[n] = dist[m] + admatrix[m][n];
	//				pre[n] = m;
	//			}
	//		}*/
	//	}
	//}
}
