#include "LoopDef.h"
namespace LoopGen
{
	cylinder::cylinder(cylinder&& rhs)
	{
		id = rhs.id;
		mesh = rhs.mesh;
		vertices = std::move(rhs.vertices);
		bounds = std::move(rhs.bounds);
		faces = std::move(rhs.faces);
		vertice_flag = std::move(rhs.vertice_flag);
		info_on_region = std::move(rhs.info_on_region);
		bound_flag = std::move(rhs.bound_flag);
		face_flag = std::move(rhs.face_flag);
		cut_v_flag = std::move(rhs.cut_v_flag);
		cut_f_flag = std::move(rhs.cut_f_flag);
		vidmap = std::move(rhs.vidmap);
		uv[0] = std::move(rhs.uv[0]);
		uv[1] = std::move(rhs.uv[1]);
	}

	void cylinder::SetBound()
	{
		std::deque<bool> bv_flag(mesh->n_vertices(), false);
		for (auto tv : vertices)
			for (auto tvv : mesh->vv_range(tv))
				if (!vertice_flag[tvv.idx()])
					bv_flag[tv.idx()] = true;
		bound_flag.resize(mesh->n_halfedges(), false);

		std::deque<bool> searched(mesh->n_vertices(), false);
		for (auto tv : vertices)
		{
			if (!bv_flag[tv.idx()] || searched[tv.idx()])
				continue;
			searched[tv.idx()] = true;
			HalfedgeHandle he = mesh->voh_begin(tv);
			while (!face_flag[mesh->face_handle(he).idx()] && face_flag[mesh->face_handle(mesh->opposite_halfedge_handle(he)).idx()])
			{
				he = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));
			}
			HalfedgeHandle he_transfer = he;
			std::vector<HalfedgeHandle> bnd;
			do
			{
				searched[mesh->to_vertex_handle(he).idx()] = true;
				bnd.push_back(he_transfer);
				he_transfer = mesh->opposite_halfedge_handle(he_transfer);
				while (!face_flag[mesh->face_handle(he_transfer).idx()] && face_flag[mesh->face_handle(mesh->opposite_halfedge_handle(he_transfer)).idx()])
				{
					he_transfer = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he_transfer));
				}
			} while (he_transfer != he);
			bounds.push_back(std::move(bnd));
		}
		if (bounds.size() != 2)
		{
			dprint("a region that does not have two boundaries");
			system("pause");
		}
	}

	OpenMesh::Vec3d cylinder::GetUGrad(FaceHandle fh)
	{
		std::vector<int> vid; vid.reserve(3);
		std::vector<double> u; u.reserve(3);
		for (auto fv : mesh->fv_range(fh))
		{
			vid.push_back(fv.idx());
			u.push_back(uv[0](vidmap[fv.idx()]));
		}
		if (cut_f_flag[fh.idx()])
		{
			for (int i = 0; i < 3; ++i)
			{
				if (cut_v_flag[vid[i]])
					u[i] += 1.0;
			}
		}
		OpenMesh::Vec3d grad(0.0, 0.0, 0.0);
		for (int i = 0; i < 3; ++i)
		{
			grad += u[i] * mesh->normal(fh).cross(mesh->calc_edge_vector(mesh->find_halfedge(mesh->vertex_handle(vid[(i + 1) % 3]), mesh->vertex_handle(vid[(i + 2) % 3]))));
		}
		return grad * 0.5 / mesh->calc_face_area(fh);
	}

	OpenMesh::Vec3d cylinder::GetVGrad(FaceHandle fh)
	{
		std::vector<int> vid; vid.reserve(3);
		std::vector<double> v; v.reserve(3);
		for (auto fv : mesh->fv_range(fh))
		{
			vid.push_back(fv.idx());
			v.push_back(uv[1](vidmap[fv.idx()]));
		}
		OpenMesh::Vec3d grad(0.0, 0.0, 0.0);
		for (int i = 0; i < 3; ++i)
		{
			grad += v[i] * mesh->normal(fh).cross(mesh->calc_edge_vector(mesh->find_halfedge(mesh->vertex_handle(vid[(i + 1) % 3]), mesh->vertex_handle(vid[(i + 2) % 3]))));
		}
		return grad * 0.5 / mesh->calc_face_area(fh);
	}
}