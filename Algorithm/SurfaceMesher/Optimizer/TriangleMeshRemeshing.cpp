#include "TriangleMeshRemeshing.h"

namespace CADMesher
{
	void TriangleMeshRemeshing::run()
	{
		for (auto tv : mesh->vertices())
			mesh->data(tv).set_targetlength(expected_length);

		printMeshQuality(*mesh);
		for (int i = 0; i < 30; i++)
		{
			dprint("\niteration times:", i + 1);
			split();
			collapse();
			equalize_valence();
			if (i > 10 && meshMinAngle(*mesh) < lowerAngleBound)
			{
				adjustTargetLength();
				//processAngle();
			}
			tangential_relaxation();
			dprint("mesh vertices number:", mesh->n_vertices());
			printMeshQuality(*mesh);
		}
		dprint("Remeshing Done!\n");

		printMeshQuality(*mesh);
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
		if (polymeshInput)
		{
			assembleMesh();
		}
#endif
	}

	void TriangleMeshRemeshing::split()
	{
		for (auto te : mesh->edges()) {
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
			if (polymeshInput && te.v0().idx() < boundaryNum && te.v1().idx() < boundaryNum)
			{
				continue;
			}
#endif
			//dprint(te.idx());
			SOV vert[2] = { te.v0(), te.v1() };
			double t0 = mesh->data(vert[0]).get_targetlength();
			double t1 = mesh->data(vert[1]).get_targetlength();
			if (mesh->calc_edge_length(te) < 1.333 * std::min(t0, t1) || te.is_boundary())
				continue;
			OV newvert = mesh->add_vertex(mesh->calc_edge_midpoint(te));
			mesh->data(newvert).set_targetlength(0.5*(t0 + t1));
			//bool flag = mesh->data(te).get_edgeflag();
			bool flag1 = mesh->data(te).flag1;
			bool flag2 = mesh->data(te).flag2;
			mesh->split_edge(te, newvert);

			if (flag1)
			{
				mesh->data(newvert).set_vertflag(true);
				mesh->data(mesh->edge_handle(mesh->find_halfedge(vert[0], newvert))).flag1 = true;
				mesh->data(mesh->edge_handle(mesh->find_halfedge(vert[1], newvert))).flag1 = true;
			}
			else if (flag2)
			{
				mesh->data(newvert).set_vertflag(true);
				mesh->data(mesh->edge_handle(mesh->find_halfedge(vert[0], newvert))).flag2 = true;
				mesh->data(mesh->edge_handle(mesh->find_halfedge(vert[1], newvert))).flag2 = true;
			}
			else
			{
				mesh->set_point(newvert, aabbtree->closest_point(GravityPos(newvert)));
			}
		}
		mesh->garbage_collection();
		//mesh->update_normals();
		dprint("split done");
	}

	void TriangleMeshRemeshing::collapse()
	{
		for (auto the : mesh->halfedges()) {
			if (!mesh->is_collapse_ok(the))
			{
				continue;
			}
#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
			if (polymeshInput && the.from().idx() < boundaryNum)
			{
				continue;
			}
#endif
			OV fromvert = mesh->from_vertex_handle(the);
			OV tovert = mesh->to_vertex_handle(the);
			double t0 = mesh->data(fromvert).get_targetlength();
			double t1 = mesh->data(tovert).get_targetlength();
			double x = mesh->calc_edge_length(the);
			double min_of_t0_t1 = std::min(t0, t1);

			if (mesh->data(fromvert).get_vertflag()) 
				continue;
			if (x >= 0.8 * min_of_t0_t1)
				continue;

#if 0
			//特征线上移动点
			if (mesh->data(the.edge()).get_edgeflag()) {
				if (x >= 0.666 * min_of_t0_t1) continue;
				int count = 0;
				OV v = fromvert;
				for (auto tvoh : mesh->voh_range(fromvert)) {
					if (!mesh->data(tvoh.edge()).get_edgeflag()) continue;
					if (tvoh.to() == tovert) continue;
					++count;
					v = tvoh.to();
				}
				if (count == 1) {
					double y = (mesh->point(v) - mesh->point(tovert)).norm();
					double z = (mesh->point(v) - mesh->point(fromvert)).norm();
					double a = sqrt((x + y + z)*(x + y - z)*(x + z - y)*(y + z - x)) * 0.25;
					//if (400 * a > expected_length*y) continue;//要保证由x,y,z构成的三角形在y上的高，小于err=expected_length*0.05
					//这里为了进一步固定特征，要求更严格
					if (z > x)
					{
						mesh->set_point(fromvert, aabbtree->closest_point((mesh->point(v) + mesh->point(tovert))*0.5));
					}
					else {
						mesh->collapse(the);
						mesh->data(mesh->edge_handle(mesh->find_halfedge(v, tovert))).set_edgeflag(true);
					}
				}
				continue;
			}
			//if (mesh->data(mesh->edge_handle(the)).get_edgeflag()) continue;
			//移除短边
			if (mesh->data(fromvert).get_vertflag()) {
				if (mesh->data(tovert).get_vertflag() && x < min_of_t0_t1*0.1)
				{
					mesh->collapse(the);
				}
				continue;
			}
#else
			if (mesh->data(the.edge()).flag1 || mesh->data(the.edge()).flag2)
				continue;
			if (mesh->is_boundary(the))
				continue;
#endif
			O3d pos = mesh->data(tovert).get_vertflag() ? mesh->point(tovert) : mesh->calc_centroid(the);

			//stop if collapsing results in long edges
			for (OV thev : mesh->vv_range(fromvert))
				if ((pos - mesh->point(thev)).norm() > 1.33 * min_of_t0_t1)
					goto goto20210523;
			for (OV thev : mesh->vv_range(tovert))
				if ((pos - mesh->point(thev)).norm() > 1.33 * min_of_t0_t1)
					goto goto20210523;

			//stop if collapsing results in small angles
			//auto t_f = mesh->opposite_halfedge_handle(the);
			//auto he = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(the)));
			auto t_f = the.opp();
			auto he = the.prev().opp().prev();
			while (he != t_f) {
				/*O3d p1 = mesh->point(mesh->from_vertex_handle(he));
				O3d p2 = mesh->point(mesh->opposite_vh(he));*/
				O3d p1 = mesh->point(he.from());
				O3d p2 = mesh->point(mesh->opposite_vh(he));
				O3d p01 = (p1 - pos).normalize();
				O3d p12 = (p2 - p1).normalize();
				O3d p20 = (pos - p2).normalize();
				if (acos(-p01.dot(p12)) < 0.1) goto goto20210523;
				if (acos(-p12.dot(p20)) < 0.1) goto goto20210523;
				if (acos(-p20.dot(p01)) < 0.1) goto goto20210523;
				//he = mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(he));
				he = he.opp().prev();
			}

			//stop if collapsing results in selfintersection
			//to be completed

			mesh->set_point(tovert, pos);
			mesh->collapse(the);
		goto20210523:;
		}
		mesh->garbage_collection();
		dprint("collapse done");
	}

	void TriangleMeshRemeshing::equalize_valence()
	{
		auto va = [&](OV v) {return mesh->valence(v) + (mesh->is_boundary(v) ? 2 : 0); };

		double PiDividedByThree = PI * 0.3333;
		double temp = PiDividedByThree - lowerAngleBound; temp *= temp;
		auto penaltyFunction = [&](double a)
		{
			double t = a - PiDividedByThree; t *= t;
			return t > temp ? t * t : t;
		};
		auto opt = [&](OpenMesh::SmartHalfedgeHandle &he) {
			return penaltyFunction(mesh->calc_sector_angle(he)) + penaltyFunction(mesh->calc_sector_angle(he.prev()))
				+ penaltyFunction(mesh->calc_sector_angle(he.next()));
		};
		for (auto te : mesh->edges()) {
			if (!mesh->is_flip_ok(te) || mesh->data(te).flag1 || mesh->data(te).flag2)
				continue;
			auto h0 = te.h0();
			auto h1 = te.h1();
			int v0 = va(te.v0());
			int v1 = va(te.v1());
			int u0 = va(mesh->opposite_vh(h0));
			int u1 = va(mesh->opposite_vh(h1));
			if (fabs(v0 - 6) + fabs(v1 - 6) + fabs(u0 - 6) + fabs(u1 - 6) <= fabs(v0 - 7) + fabs(v1 - 7) + fabs(u0 - 5) + fabs(u1 - 5)) continue;

			//在一定程度上可防止翻折
			double alpha0, alpha1;
			alpha0 = acos(-mesh->calc_edge_vector(h0).dot(mesh->calc_edge_vector(/*mesh->next_halfedge_handle(h0)*/h0.next())));
			alpha1 = acos(-mesh->calc_edge_vector(h1).dot(mesh->calc_edge_vector(/*mesh->prev_halfedge_handle(h1)*/h1.prev())));
			if (alpha0 + alpha1 > PI) continue;
			alpha0 = acos(-mesh->calc_edge_vector(h0).dot(mesh->calc_edge_vector(/*mesh->prev_halfedge_handle(h0)*/h0.prev())));
			alpha1 = acos(-mesh->calc_edge_vector(h1).dot(mesh->calc_edge_vector(/*mesh->next_halfedge_handle(h1)*/h1.next())));
			if (alpha0 + alpha1 > PI) continue;

			////检查二面角
			//auto n0 = mesh->calc_face_normal(/*mesh->face_handle(h0)*/h0.face());
			//auto n1 = mesh->calc_face_normal(/*mesh->face_handle(h1)*/h1.face());
			//if (n0.dot(n1) < 0.8)
			//{
			//	if (mesh->data(te.v0()).get_vertflag() && mesh->data(te.v1()).get_vertflag())
			//	{
			//		mesh->data(te).set_edgeflag(true);
			//		continue;
			//	}
			//}
			////防止出现狭长三角形
			//auto V0 = mesh->point(te->v0());
			//auto V1 = mesh->point(te->v1());
			//auto U0 = mesh->point(mesh->opposite_vh(h0));
			//auto U1 = mesh->point(mesh->opposite_vh(h1));
			//if (((U0 - V1).norm() + (U1 - V1).norm()) / (U0 - U1).norm() < 1.1) continue;
			//if (((U0 - V0).norm() + (U1 - V0).norm()) / (U0 - U1).norm() < 1.1) continue;

			//假设flip，检查局部网格角度是否被优化
			double opt_before = opt(h0) + opt(h1);
			mesh->flip(te);
			//若局部网格角度未被优化，则再次flip回到初始状态
			if (opt_before < opt(te.h0()) + opt(te.h1()))
				mesh->flip(te);
		}
		mesh->garbage_collection();
		dprint("equalize done");
	}

	void TriangleMeshRemeshing::adjustTargetLength()
	{
		double maxL = 0;
		double minL = expected_length;
		double sum = 0;
		double threshold = 0.1 * expected_length;

		for (auto tv : mesh->vertices())
		{
			sum = 0;
			for (auto tve : mesh->ve_range(tv))
			{
				sum += mesh->calc_edge_length(tve);
			}
			mesh->data(tv).set_targetlength(std::min(expected_length, std::max(threshold, 1.2*sum / mesh->valence(tv))));
		}
	}

	void TriangleMeshRemeshing::processAngle()
	{
		int id = 0;
		double threshold = 0.1*expected_length;
		double minAngle = 4.0;
		for (auto tf : mesh->faces())
		{
			if (!tf.is_valid()) continue;
			id = 0;
			minAngle = 4.0;
			int i = 0;
			for (auto tfh : mesh->fh_range(tf))
			{
				double angle = mesh->calc_sector_angle(tfh);
				if (angle < minAngle)
				{
					minAngle = angle;
					id = i;
				}
				++i;
			}
			if (minAngle < lowerAngleBound)
			{
				auto h_iter = mesh->fh_begin(tf);
				for (i = 0; i < id; ++i)
				{
					++h_iter;
				}
				auto th = mesh->prev_halfedge_handle(*h_iter);
				auto te = mesh->edge_handle(th);
				if (mesh->is_boundary(th.from()))
				{
					continue;
				}
				if (mesh->is_collapse_ok(th) && mesh->calc_edge_length(te) < threshold)
				{
					if (mesh->data(th.from()).get_vertflag())
					{
						mesh->data(th.to()).set_vertflag(true);
					}
					if (mesh->data(th.prev().edge()).get_edgeflag())
					{
						mesh->data(h_iter->edge()).set_edgeflag(true);
					}
					mesh->collapse(th);
				}
				else
				{
					if (mesh->calc_sector_angle(th) < mesh->calc_sector_angle(th.prev()))
					{
						th = th.prev();
					}
					auto flagvert = th.to();
					auto ph = th.prev();
					auto ne = th.next().edge();
					te = th.edge();
					auto pe = ph.edge();

					if (!mesh->is_flip_ok(pe))
						continue;
					if (mesh->data(pe).get_edgeflag())
					{
						mesh->data(flagvert).set_vertflag(true);
						mesh->data(te).set_edgeflag(true);
						mesh->data(ne).set_edgeflag(true);
						mesh->data(pe).set_edgeflag(false);
						/*if (!mesh->data(flagvert).get_vertflag())
						{
							mesh->set_point(flagvert, mesh->calc_edge_midpoint(pe));
						}*/
					}
					mesh->flip(pe);
				}
			}
		}
		mesh->garbage_collection();
	}

	void TriangleMeshRemeshing::tangential_relaxation()
	{
		for (auto tv : mesh->vertices()) {
			if (mesh->data(tv).get_vertflag() || mesh->is_boundary(tv))
				continue;
			//mesh->set_point(tv, GravityPos(tv));
			mesh->set_point(tv, aabbtree->closest_point(GravityPos(tv)));
		}
		dprint("project done");
	}

	O3d TriangleMeshRemeshing::GravityPos(const OV &v)
	{
		if (mesh->is_boundary(v))
		{
			return mesh->point(v);
		}
		O3d sum(0, 0, 0);
		double area = 0;
		for (auto tvoh : mesh->voh_range(v))
		{
			double a = mesh->calc_face_area(mesh->face_handle(tvoh)) + mesh->calc_face_area(mesh->opposite_face_handle(tvoh));
			area += a;
			sum += a*mesh->point(mesh->to_vertex_handle(tvoh));
		}
		if (area < epsilonerror)
		{
			for (auto tvoh : mesh->voh_range(v))
			{
				sum += mesh->point(mesh->to_vertex_handle(tvoh));
			}
			return sum / mesh->valence(v);
		}
		return sum / area;
	}

#ifdef OPENMESH_POLY_MESH_ARRAY_KERNEL_HH
	TriangleMeshRemeshing::TriangleMeshRemeshing(PolyMesh *mesh_, double target_length)
		:expected_length(target_length), polymeshInput(true)
	{
		initMeshStatusAndNormal(*mesh_);
		polymesh = new Mesh();
		polymesh->reserve(mesh_->n_vertices(), mesh_->n_edges(), mesh_->n_faces());

		enum vertexType { tri, poly, mixed };
		OpenMesh::VPropHandleT<int> idMap;
		OpenMesh::VPropHandleT<vertexType> typeMap;
		mesh_->add_property(idMap);
		mesh_->add_property(typeMap);
		int i = 0, j = 0, k = 0;
		{
			bool allTri = true, allPoly = true;
			int triangleNum = 0;
			for (auto &tv : mesh_->vertices())
			{
				allTri = allPoly = true;
				for (auto &tvf : mesh_->vf_range(tv))
				{
					if (tvf.valence() == 3)
					{
						allPoly = false;
						++triangleNum;
					}
					else
					{
						allTri = false;
					}
				}
				if (allTri)
				{
					mesh_->property(idMap, tv) = i++;
					mesh_->property(typeMap, tv) = tri;
				}
				else if (allPoly)
				{
					mesh_->property(idMap, tv) = j++;
					mesh_->property(typeMap, tv) = poly;
				}
				else
				{
					mesh_->property(idMap, tv) = k++;
					mesh_->property(typeMap, tv) = mixed;
				}
			}
			triangleNum /= 3;
			for (auto &tv : mesh_->vertices())
			{
				if (mesh_->property(typeMap, tv) != mixed)
				{
					mesh_->property(idMap, tv) += k;
				}
			}
			mesh = new TriMesh();
			mesh->reserve(i + k, i + k + triangleNum, triangleNum);//V+F-E=2-2g-b  =>  E~=V+F
		}

		for (auto &tv : mesh_->vertices())
		{
			if (mesh_->property(typeMap, tv) == mixed)
			{
				auto &av = mesh->add_vertex(mesh_->point(tv));
				auto &pav = polymesh->add_vertex(mesh_->point(tv));
				if (mesh_->data(tv).get_vertflag())
				{
					mesh->data(av).set_vertflag(true);
					polymesh->data(pav).set_vertflag(true);
				}
			}
		}
		for (auto &tv : mesh_->vertices())
		{
			switch (mesh_->property(typeMap, tv))
			{
			case tri:
			{
				auto &av = mesh->add_vertex(mesh_->point(tv));
				if (mesh_->data(tv).get_vertflag())
				{
					mesh->data(av).set_vertflag(true);
				}
				break;
			}
			case poly:
			{
				auto &av = polymesh->add_vertex(mesh_->point(tv));
				if (mesh_->data(tv).get_vertflag())
				{
					polymesh->data(av).set_vertflag(true);
				}
				break;
			}
			default:
				break;
			}
		}

		int id[4];
		int ii = 0;
		for (auto &tf : mesh_->faces())
		{
			ii = 0;
			if (tf.valence() == 3)
			{
				for (auto &tfv : mesh_->fv_range(tf))
				{
					id[ii++] = mesh_->property(idMap, tfv);
				}
				mesh->add_face(mesh->vertex_handle(id[0]), mesh->vertex_handle(id[1]), mesh->vertex_handle(id[2]));
			}
			else
			{
				for (auto &tfv : mesh_->fv_range(tf))
				{
					id[ii++] = mesh_->property(idMap, tfv);
				}
				polymesh->add_face(polymesh->vertex_handle(id[0]), polymesh->vertex_handle(id[1]),
					polymesh->vertex_handle(id[2]), polymesh->vertex_handle(id[3]));
			}
		}
		for (auto &te : mesh_->edges())
		{
			if (mesh_->data(te).get_edgeflag())
			{
				vertexType vt0 = mesh_->property(typeMap, te.v0());
				vertexType vt1 = mesh_->property(typeMap, te.v1()); 
				auto th = mesh_->find_halfedge(mesh_->vertex_handle(mesh_->property(idMap, te.v0())),
					mesh_->vertex_handle(mesh_->property(idMap, te.v0())));
				if (vt0 == tri || vt1 == tri)
				{
					mesh->data(mesh->edge_handle(mesh->find_halfedge(mesh->vertex_handle(mesh_->property(idMap, te.v0())),
						mesh->vertex_handle(mesh_->property(idMap, te.v1()))))).set_edgeflag(true);
				}
				else if (vt0 == poly || vt1 == poly)
				{
					polymesh->data(polymesh->edge_handle(polymesh->find_halfedge(polymesh->vertex_handle(mesh_->property(idMap, te.v0())),
						polymesh->vertex_handle(mesh_->property(idMap, te.v1()))))).set_edgeflag(true);
				}
				else
				{
					mesh->data(mesh->edge_handle(mesh->find_halfedge(mesh->vertex_handle(mesh_->property(idMap, te.v0())),
						mesh->vertex_handle(mesh_->property(idMap, te.v1()))))).set_edgeflag(true);
					polymesh->data(polymesh->edge_handle(polymesh->find_halfedge(polymesh->vertex_handle(mesh_->property(idMap, te.v0())),
						polymesh->vertex_handle(mesh_->property(idMap, te.v1()))))).set_edgeflag(true);
				}
			}
		}
		boundaryNum = k;
		mesh_->remove_property(idMap);
		mesh_->remove_property(typeMap);

		*mesh_ = *polymesh;
		delete polymesh;
		polymesh = mesh_;
		initMeshStatusAndNormal(*mesh);

		//*polymesh = Mesh(*mesh);

		if (expected_length <= 0)
		{
			expected_length = meshAverageLength(*mesh);
		}
		high = 1.33*expected_length;
		low = 0.8*expected_length;
		aabbtree = new ClosestPointSearch::AABBTree(*mesh);
	}

	void TriangleMeshRemeshing::assembleMesh()
	{
#if 1
		int nv = polymesh->n_vertices();
		auto vItr = mesh->vertices_begin();
		for (int i = 0; i < boundaryNum; ++i, ++vItr);
		for (; vItr != mesh->vertices_end(); ++vItr)
		{
			auto &av = polymesh->add_vertex(mesh->point(*vItr));
			if (mesh->data(*vItr).get_vertflag())
			{
				polymesh->data(av).set_vertflag(true);
			}
		}

		int id[3];
		int ii = 0;
		for (auto &tf : mesh->faces())
		{
			ii = 0;
			for (auto &tfv : mesh->fv_range(tf))
			{
				id[ii++] = tfv.idx() < boundaryNum ? tfv.idx() : (tfv.idx() + nv - boundaryNum);
			}
			polymesh->add_face(polymesh->vertex_handle(id[0]), polymesh->vertex_handle(id[1]), polymesh->vertex_handle(id[2]));
		}
		for (auto &te : mesh->edges())
		{
			if (mesh->data(te).get_edgeflag())
			{
				polymesh->data(polymesh->edge_handle(polymesh->find_halfedge(
					polymesh->vertex_handle(te.v0().idx() < boundaryNum ? te.v0().idx() : te.v0().idx() + nv),
					polymesh->vertex_handle(te.v1().idx() < boundaryNum ? te.v1().idx() : te.v1().idx() + nv)))).set_edgeflag(true);
			}
		}
#else
		*polymesh = Mesh(*mesh);
#endif
		initMeshStatusAndNormal(*polymesh);

		delete mesh;
		mesh = nullptr;
	}
#endif

	/*bool TriangleMeshRemeshing::split_one_edge(Mesh::EdgeHandle& eh, OpenMesh::Vec3d& p)
	{
		Mesh::HalfedgeHandle heh0 = mesh->halfedge_handle(eh, 0);
		Mesh::HalfedgeHandle heh1 = mesh->halfedge_handle(eh, 1);
		Mesh::VertexHandle vh0 = mesh->to_vertex_handle(heh0); OpenMesh::Vec3d p0 = mesh->point(vh0);
		Mesh::VertexHandle vh1 = mesh->to_vertex_handle(heh1); OpenMesh::Vec3d p1 = mesh->point(vh1);

		std::vector<Mesh::VertexHandle> one_face(3);
		bool flag = mesh->data(eh).get_edgeflag();
		if (mesh->is_boundary(eh))
		{
			if (Mesh::InvalidFaceHandle != mesh->face_handle(heh0))
			{
				Mesh::VertexHandle vh2 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh0));
				OpenMesh::Vec3d p2 = mesh->point(vh2);
				OpenMesh::Vec3d n = OpenMesh::cross(p1 - p2, p0 - p2).normalize();
				double a1 = OpenMesh::dot(n, OpenMesh::cross(p2 - p0, p - p0));
				double a2 = OpenMesh::dot(n, OpenMesh::cross(p1 - p2, p - p2));
				if (a1 < 1e-8 || a2 < 1e-8) return false;
				Mesh::VertexHandle vh = mesh->add_vertex(p);
				mesh->delete_edge(eh, false); mesh->garbage_collection();
				one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh->add_face(one_face);
				one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh->add_face(one_face);
			}
			else
			{
				Mesh::VertexHandle vh3 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh1));
				OpenMesh::Vec3d p3 = mesh->point(vh3);
				OpenMesh::Vec3d n = OpenMesh::cross(p0 - p3, p1 - p3).normalize();
				double a1 = OpenMesh::dot(n, OpenMesh::cross(p0 - p3, p - p3));
				double a2 = OpenMesh::dot(n, OpenMesh::cross(p3 - p1, p - p1));
				if (a1 < 1e-8 || a2 < 1e-8) return false;
				Mesh::VertexHandle vh = mesh->add_vertex(p);
				mesh->delete_edge(eh, false); mesh->garbage_collection();
				one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh->add_face(one_face);
				one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh->add_face(one_face);
			}
		}
		else
		{
			Mesh::VertexHandle vh2 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh0)); OpenMesh::Vec3d p2 = mesh->point(vh2);
			OpenMesh::Vec3d n1 = OpenMesh::cross(p1 - p2, p0 - p2).normalize();
			double a1 = OpenMesh::dot(n1, OpenMesh::cross(p2 - p0, p - p0));
			double a2 = OpenMesh::dot(n1, OpenMesh::cross(p1 - p2, p - p2));
			Mesh::VertexHandle vh3 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh1)); OpenMesh::Vec3d p3 = mesh->point(vh3);
			OpenMesh::Vec3d n2 = OpenMesh::cross(p0 - p3, p1 - p3).normalize();
			double a3 = OpenMesh::dot(n2, OpenMesh::cross(p0 - p3, p - p3));
			double a4 = OpenMesh::dot(n2, OpenMesh::cross(p3 - p1, p - p1));
			if (a1 < 1e-8 || a2 < 1e-8 || a3 < 1e-8 || a4 < 1e-8) return false;
			Mesh::VertexHandle vh = mesh->add_vertex(p);
			mesh->delete_edge(eh, false); mesh->garbage_collection();
			one_face[0] = vh0; one_face[1] = vh2; one_face[2] = vh; mesh->add_face(one_face);
			one_face[0] = vh2; one_face[1] = vh1; one_face[2] = vh; mesh->add_face(one_face);
			one_face[0] = vh3; one_face[1] = vh0; one_face[2] = vh; mesh->add_face(one_face);
			one_face[0] = vh1; one_face[1] = vh3; one_face[2] = vh; mesh->add_face(one_face);
		}
		return true;
	}*/
}