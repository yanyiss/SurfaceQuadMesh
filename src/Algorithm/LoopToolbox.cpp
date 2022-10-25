#include "LoopToolbox.h"
#include <Eigen\Cholesky>
#include <omp.h>
namespace LoopGen
{
	bool LoopGen::FieldAligned_PlanarLoop(VertexHandle v, std::vector<VertexHandle>& loop, int shift)
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

			//loop.push_back(vert.id); loop.push_back(mesh->from_vertex_handle(prev[vert.id]).idx());
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
					//double dot_ = (position.col(toid) - position.col(fromid)).normalized().dot(plane_normal); dot_ *= dot_;
					//w = sqrt(w + dot_);
					w = sqrt(w * 900 + 1 - w);
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
#if 0
		cf = new crossField(mesh);
#else
		cf = new crossField("..//resource//field//vase.field");
#endif
		dprint("Initialize Field Done!");
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

			/*if (eitr->idx() == 6097)
			{
				dprint("eitr0");
				dprint(fv.dot(ev.normalized()), gv.dot(ev.normalized()), crossfield.col(4 * gid).dot(ev.normalized()));
				dprint(mesh->from_vertex_handle(h0).idx(), mesh->to_vertex_handle(h0).idx());
				dprint(matching[h0.idx()]);
				int p = 0;
			}*/

			if (s < c)
				c = DBL_MAX;
			else
				s = DBL_MAX;
			if (arc >= 0 && arc < halfPI)
				w0 << s, DBL_MAX, DBL_MAX, c;
			else if (arc >= halfPI && arc < PI)
				w0 << DBL_MAX, DBL_MAX, s, c;
			else if (arc >= _PI && arc < _halfPI)
				w0 << DBL_MAX, c, s, DBL_MAX;
			else
				w0 << s, c, DBL_MAX, DBL_MAX;

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
		eov.resize(mesh->n_vertices());
		auto start = clock();
#pragma omp parallel for
		for (int i = 0; i < mesh->n_vertices(); ++i)
		{
			std::vector<VertexHandle> loop;
			double plane[4];
			Eigen::VectorXd xyz[3];
			if (!FieldAligned_PlanarLoop(mesh->vertex_handle(i), loop, 0))
			{
				eov[i] = DBL_MAX;
			}
			else
			{
				GetPositionFromLoop(loop, xyz);
				eov[i] = EvaluatePlanarity(xyz, plane);
			}
			if (!FieldAligned_PlanarLoop(mesh->vertex_handle(i), loop, 1))
			{
				eov[i] = std::min(eov[i], DBL_MAX);
			}
			else
			{
				GetPositionFromLoop(loop, xyz);
				eov[i] = std::min(EvaluatePlanarity(xyz, plane), eov[i]);
			}
			//dprint("index:", i, eov[i]);
		}
		dprint("time:", clock() - start);
		std::ofstream file_writer;
		file_writer.open("..//resource//energy//vase.energy");
		if (file_writer.fail()) {
			std::cout << "fail to open\n";
		}
		for (auto e : eov)
			file_writer << e << "\n";
		file_writer.close();
	}

	bool LoopGen::RefineLoop(std::vector<VertexHandle>& loop, PlaneLoop& planar_loop)
	{
		Eigen::VectorXd xyz[3];
		GetPositionFromLoop(loop, xyz);
		double plane[4];
		if (EvaluatePlanarity(xyz, plane) > 0.01)
			//return false;
			int p = 0;
		auto dis = [&](Vec3d& pos)
		{
			return plane[0] * pos[0] + plane[1] * pos[1] + plane[2] * pos[2] + plane[3];
		};

		//从起始点出发，不断加入端点落在平面两侧的边
		auto hitr = mesh->voh_begin(loop[0]);
		auto s0 = dis(mesh->point(mesh->to_vertex_handle(mesh->next_halfedge_handle(hitr.handle()))));
		double distance[2];
		PointOnHalfedge poh[2]; int id = 0;
		for (; hitr != mesh->voh_end(loop[0]); ++hitr)
		{
			auto s1 = dis(mesh->point(mesh->to_vertex_handle(hitr.handle())));
			if (s0 * s1 < 0)
			{
				/*h[id] = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(hitr.handle()));
				++id;*/
				if (s0 > 0)
				{
					poh[0].h = mesh->next_halfedge_handle(hitr.handle());
					poh[0].c = s0 / (s0 - s1);
					distance[0] = s0; distance[1] = s1;
				}
				else
				{
					poh[1].h = mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(hitr.handle()));
					poh[1].c = s1 / (s1 - s0);
				}
				++id;
			}
			s0 = s1;
		}
		if (id != 2)
		{
			dprint("error out");
			system("pause");
		}

		planar_loop.clear();
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
		}
		return true;
#if 0
		//首先向平面投影，再向网格投影
		Eigen::VectorXd v(xyz[0](0), xyz[1](0), xyz[2](0));
		Eigen::VectorXd n(plane[0], plane[1], plane[2]);

		int nv = loop.size() - 1;
		planar_loop.resize(3, nv);
		for (int i = 0; i < nv; ++i)
		{
			auto& pos = mesh->point(loop[i]);
			Eigen::VectorXd p(pos[0], pos[1], pos[2]);
			auto plane_pos = n.dot(p - v) * n + p;
			auto proj_pos = aabbtree->closest_point_and_face_handle(Vec3d(plane_pos(0), plane_pos(1), plane_pos(2)));
			planar_loop.col(i) << proj_pos.first[0], proj_pos.first[1], proj_pos.first[2];
		}
#endif
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
		LeastSquarePlane(xyz, plane);
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
		bounding[0] = DBL_MAX; bounding[1] = -DBL_MAX;
		bounding[2] = DBL_MAX; bounding[3] = -DBL_MAX;
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