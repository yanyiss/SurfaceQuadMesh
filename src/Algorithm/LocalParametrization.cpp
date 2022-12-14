#include "LocalParametrization.h"
namespace LoopGen
{
	LocalParametrization::LocalParametrization(Mesh& mesh_, VertexHandle v)
	{
		mesh = &mesh_;
		int nf = mesh->n_faces();
		int nv = mesh->n_vertices();
		region_f_flag.resize(nf, false);
		region_v_flag.resize(nv, false);
		region_vertex.push_back(v);
		region_v_flag[v.idx()] = true;
		uv[0].resize(1); uv[0].setZero();
		uv[1].resize(1); uv[1].setZero();
		all_pl.resize(nv);
		vidmap.resize(nv); vidmap[v.idx()] = 0;
		x_axis.resize(3, nf);
		y_axis.resize(3, nf);
		grow_dir.resize(nv, -1);
	}

	void LocalParametrization::run(const Eigen::Matrix3Xd& normal)
	{
		int nf = mesh->n_faces();

		//标记与cut相关的顶点和面，计算装配矩阵需要的数据
		int nv = mesh->n_vertices();
		cutv_flag.resize(nv, false);
		cutf_flag.resize(nf, false);
		std::vector<std::map<VertexHandle, std::pair<bool, Eigen::Vector3d>>> info(nf);
		{
			for (auto c : cut)
				cutv_flag[c.idx()] = true;

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
			{
				vidmap[vv.idx()] = count++;
			}
		}

		int new_vertex_size = new_vertex.size();
		int new_face_size = new_face.size();
		std::vector<Eigen::Triplet<double>> triple;
		std::vector<double> w(nv, 0);
		std::vector<double> size_ratio(nf, 1.0);
		uv[0].conservativeResize(region_vertex_size + new_vertex_size); uv[0].tail(new_vertex_size).setZero();
		uv[1].conservativeResize(region_vertex_size + new_vertex_size); uv[1].tail(new_vertex_size).setZero();
		//Eigen::Vector3d right[2];
		//right[0].resize(new_vertex_size + new_face_size); right[0].setZero();
		//right[1].resize(new_vertex_size + new_face_size); right[1].setZero();

		//const auto& crossfield = cf->getCrossField();
		//for (auto v : new_vertex)
			//dprint(v.idx(), vidmap[v.idx()]);
		for (auto v : new_vertex)
		{
			int vid = v.idx();
			int vm = vidmap[vid];
			//dprint(vid, vm);
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
					if (f_info.second.first)
						uv[0](vm) -= dot_;
					//right[0](vm) -= dot_;
					if (!new_v_flag[fvid])
					{
						uv[0](vm) -= dot_ * GetU(fvid);
						uv[1](vm) -= dot_ * GetV(fvid);
						//right[0](vm) -= dot_ * GetU(fvid);
						//right[1](vm) -= dot_ * GetV(fvid);
					}
					else
						w[fvid] += dot_;
				}
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
}