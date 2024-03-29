#include "LocalParametrization.h"
namespace LoopGen
{
	LocalParametrization::LocalParametrization(VertexLayer* vl_, cylinder &cy_, spread_info &sp_)
	{
		cy = &cy_;
		sp = &sp_;
		int nfl = sp->m4->facelayers.size();
		int nvl = sp->m4->verticelayers.size();
		cy->face_flag.resize(nfl, false);
		cy->vertice_flag.resize(nvl, false);
		cy->vertices.push_back(vl_);
		cy->vertice_flag[vl_->id] = true;
		cy->uv[0].resize(1); cy->uv[0].setZero();
		cy->uv[1].resize(1); cy->uv[1].setZero();
		cy->vidmap.resize(nvl, -1); cy->vidmap[vl_->id] = 0;

		sp->all_pl.resize(nvl);
		sp->x_axis.resize(3, sp->m4->mesh->n_faces()); sp->x_axis.setZero();
		sp->y_axis.resize(3, sp->m4->mesh->n_faces()); sp->y_axis.setZero();
		sp->grow_dir.resize(nvl, -1);
	}

	void LocalParametrization::run(const Eigen::Matrix3Xd& normal)
	{
		//int nf = mesh->n_faces();
		M4* m4 = sp->m4;
		auto &new_vertex = sp->new_vertex;
		auto &new_v_flag = sp->new_v_flag;
		auto &region_vertex = cy->vertices;
		auto &cut = cy->cut;
		auto &new_face = sp->new_face;
		auto &new_f_flag = sp->new_f_flag;
		auto &vidmap = cy->vidmap;
		auto &x_axis = sp->x_axis;
		auto &y_axis = sp->y_axis;
		Eigen::VectorXd* uv = cy->uv;
		int nfl = m4->facelayers.size();

		//标记与cut相关的顶点和面，计算装配矩阵需要的数据
		int nvl = m4->verticelayers.size();
		cutv_flag.resize(nvl, false);
		cutf_flag.resize(nfl, false);
		std::vector<std::map<VertexLayer*, std::pair<bool, Eigen::Vector3d>>> info(nfl);
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
		std::vector<double> w(nvl, 0);
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
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
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

	void LocalParametrization::modify_cut()
	{
		auto &cut = cy->cut;
		auto &region_f_flag = cy->face_flag;
		if (cut.size() < 3)
			return;

		int first = 0;
		while (true)
		{
			HalfedgeLayer* hl_begin = cut[first]->hl;
			HalfedgeLayer* hl_transfer = hl_begin;
			do
			{
				if (region_f_flag[hl_transfer->left])
					goto goto0;
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			++first;
		}
	goto0:;

		int last = cut.size() - 1;
		while (true)
		{
			HalfedgeLayer* hl_begin = cut[last]->hl;
			HalfedgeLayer* hl_transfer = hl_begin;
			do
			{
				if (region_f_flag[hl_transfer->left])
					goto goto1;
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			--last;
		}
	goto1:;

		for (int i = 0; i < first; ++i)
			cutv_flag[cut[i]->id] = false;
		for (int i = cut.size() - 1; i > last; --i)
			cutv_flag[cut[i]->id] = false;

		std::vector<VertexLayer*> cut_temp;
		for (int i = first; i <= last; ++i)
			cut_temp.push_back(cut[i]);
		cut = std::move(cut_temp);
	}
}