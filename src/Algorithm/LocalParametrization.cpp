#include "LocalParametrization.h"
namespace LoopGen
{
	LocalParametrization::LocalParametrization(M4& m4_, VertexLayer* vl_)
	{
		m4 = &m4_;
		int nf = m4->facelayers.size();
		int nv = m4->verticelayers.size();
		region_f_flag.resize(nf, false);
		region_v_flag.resize(nv, false);
		region_vertex.push_back(vl_);
		region_v_flag[vl_->id] = true;
		uv[0].resize(1); uv[0].setZero();
		uv[1].resize(1); uv[1].setZero();
		all_pl.resize(nv);
		vidmap.resize(nv); vidmap[vl_->id] = 0;
		x_axis.resize(3, m4->mesh->n_faces());
		y_axis.resize(3, m4->mesh->n_faces());
		grow_dir.resize(nv, -1);

	}

	void LocalParametrization::run(const Eigen::Matrix3Xd& normal)
	{
		//int nf = mesh->n_faces();
		int nf = m4->facelayers.size();

		//标记与cut相关的顶点和面，计算装配矩阵需要的数据
		int nv = m4->verticelayers.size();
		cutv_flag.resize(nv, false);
		cutf_flag.resize(nf, false);
		std::vector<std::map<VertexLayer*, std::pair<bool, Eigen::Vector3d>>> info(nf);
		{
			for (auto c : cut)
				cutv_flag[c->id] = true;
			HalfedgeLayer* hl = m4->find_halfedge_layer(cut[0], cut[1]);
			FaceLayer* fl = &m4->facelayers[hl->left];
			int flid = fl->id;
			int fid = fl->f.idx();
			Mesh* mesh = m4->mesh;
			auto calc_vector = [&](bool flag)
			{
				double inv_area = 1.0 / (2 * mesh->calc_face_area(fl->f));
				auto vi = hl->to;
				auto ev = mesh->calc_edge_vector(mesh->prev_halfedge_handle(hl->h));
				info[flid].insert(std::make_pair(&m4->verticelayers[vi], 
					std::make_pair(flag && cutv_flag[vi], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
				vi = hl->next->to;
				ev = mesh->calc_edge_vector(hl->h);
				info[flid].insert(std::make_pair(&m4->verticelayers[vi], 
					std::make_pair(flag && cutv_flag[vi], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
				vi = hl->from;
				ev = mesh->calc_edge_vector(mesh->next_halfedge_handle(hl->h));
				info[flid].insert(std::make_pair(&m4->verticelayers[vi], 
					std::make_pair(flag && cutv_flag[vi], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
			};

			while (new_f_flag[flid])
			{
				cutf_flag[flid] = true;
				calc_vector(true);
				//he = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));
				//f = mesh->face_handle(he);
				//fid = f.idx();
				hl = hl->prev->oppo;
				fl = &m4->facelayers[hl->left];
				flid = fl->id;
				fid = fl->f.idx();
			}
			for (int i = 1; i < cut.size() - 1; ++i)
			{
				//he = mesh->find_halfedge(cut[i], cut[i + 1]);
				//f = mesh->face_handle(he);
				//fid = f.idx();
				hl = m4->find_halfedge_layer(cut[i], cut[i + 1]);
				fl = &m4->facelayers[hl->left];
				flid = fl->id;
				fid = fl->f.idx();
				do
				{
					if (info[flid].empty() && new_f_flag[flid])
					{
						cutf_flag[flid] = true;
						calc_vector(true);
					}
					//he = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(he));
					//f = mesh->face_handle(he);
					//fid = f.idx();
					hl = hl->prev->oppo;
					fl = &m4->facelayers[hl->left];
					flid = fl->id;
					fid = fl->f.idx();
				} while (hl->to != cut[i - 1]->id);
			}
			//he = mesh->next_halfedge_handle(mesh->find_halfedge(cut[cut.size() - 2], cut.back()));
			//f = mesh->face_handle(he);
			//fid = f.idx();
			hl = m4->find_halfedge_layer(cut[cut.size() - 2], cut.back())->next;
			fl = &m4->facelayers[hl->left];
			flid = fl->id;
			fid = fl->f.idx();
			while (new_f_flag[flid])
			{
				if (info[flid].empty())
				{
					cutf_flag[flid] = true;
					calc_vector(true);
				}
				//he = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(he));
				//f = mesh->face_handle(he);
				//fid = f.idx();
				hl = hl->oppo->next;
				fl = &m4->facelayers[hl->left];
				flid = fl->id;
				fid = fl->f.idx();
			}

			for (auto fa : new_face)
			{
				if (cutf_flag[fa->id])
					continue;
				//he = mesh->fh_begin(fa).handle();
				//f = fa;
				//fid = f.idx();
				hl = fa->hl;
				fl = fa;
				flid = fl->id;
				fid = fl->f.idx();
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
				vidmap[vv->id] = count++;
			}
		}

		int new_vertex_size = new_vertex.size();
		int new_face_size = new_face.size();
		std::vector<Eigen::Triplet<double>> triple;
		std::vector<double> w(nv, 0);
		//std::vector<double> size_ratio(nf, 1.0);
		uv[0].conservativeResize(region_vertex_size + new_vertex_size); uv[0].tail(new_vertex_size).setZero();
		uv[1].conservativeResize(region_vertex_size + new_vertex_size); uv[1].tail(new_vertex_size).setZero();

		for (auto vl : new_vertex)
		{
			int vlid = vl->id;
			int vm = vidmap[vlid];
			//dprint(vid, vm);
			HalfedgeLayer* hl_begin = vl->hl;
			HalfedgeLayer* hl_transfer = hl_begin;
			do
			{
				int vfl_id = hl_transfer->left;
				if (new_f_flag[vfl_id])
				{
					Eigen::Vector3d& R = info[vfl_id][vl].second;
					for (const auto& f_info : info[vfl_id])
					{
						double dot_ = R.dot(f_info.second.second);
						int fvlid = f_info.first->id;
						if (f_info.second.first)
							uv[0](vm) -= dot_;
						//right[0](vm) -= dot_;
						if (!new_v_flag[fvlid])
						{
							uv[0](vm) -= dot_ * GetU(fvlid);
							uv[1](vm) -= dot_ * GetV(fvlid);
							//right[0](vm) -= dot_ * GetU(fvid);
							//right[1](vm) -= dot_ * GetV(fvid);
						}
						else
							w[fvlid] += dot_;
					}
					int vf_id = vfl_id / 4;
					uv[0](vm) += x_axis.col(vf_id).dot(R);// * size_ratio[vf_id];
					uv[1](vm) += y_axis.col(vf_id).dot(R);// *size_ratio[vf_id];
					//uv[0](vm) += crossfield.col(vf_id).dot(R) * size_ratio[vf_id];
					//uv[1](vm) += crossfield.col(vf_id % 4 == 3 ? vf_id - 3 : vf_id + 1).dot(R) * size_ratio[vf_id];
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			//for (auto vf = mesh->vf_begin(v); vf != mesh->vf_end(v); ++vf)
			//{
			//	int vf_id = vf->idx();
			//	if (!new_f_flag[vf_id])
			//		continue;
			//	Eigen::Vector3d& R = info[vf_id][v].second;
			//	for (const auto& f_info : info[vf_id])
			//	{
			//		double dot_ = R.dot(f_info.second.second);
			//		int fvid = f_info.first.idx();
			//		if (f_info.second.first)
			//			uv[0](vm) -= dot_;
			//		//right[0](vm) -= dot_;
			//		if (!new_v_flag[fvid])
			//		{
			//			uv[0](vm) -= dot_ * GetU(fvid);
			//			uv[1](vm) -= dot_ * GetV(fvid);
			//			//right[0](vm) -= dot_ * GetU(fvid);
			//			//right[1](vm) -= dot_ * GetV(fvid);
			//		}
			//		else
			//			w[fvid] += dot_;
			//	}
			//	uv[0](vm) += x_axis.col(vf_id).dot(R) * size_ratio[vf_id];
			//	uv[1](vm) += y_axis.col(vf_id).dot(R) * size_ratio[vf_id];
			//}
			triple.emplace_back(vm - region_vertex_size, vm - region_vertex_size, w[vlid]);
			w[vlid] = 0;

			hl_transfer = hl_begin;
			do
			{
				int vvid = hl_transfer->to;
				if (new_v_flag[vvid])
				{
					triple.emplace_back(vm - region_vertex_size, vidmap[vvid] - region_vertex_size, w[vvid]);
					w[vvid] = 0;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			/*for (auto vv = mesh->vv_begin(v); vv != mesh->vv_end(v); ++vv)
			{
				int vvid = vv->idx();
				if (!new_v_flag[vvid])
					continue;
				triple.emplace_back(vm - region_vertex_size, vidmap[vvid] - region_vertex_size, w[vvid]);
				w[vvid] = 0;
			}*/
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