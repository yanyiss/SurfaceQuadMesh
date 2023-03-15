#include "LoopDef.h"
namespace LoopGen
{
	void cylinder::set_bound()
	{
		//return;
		for (int i = 0; i < 2; ++i)
		{
			std::vector<HalfedgeLayer*> one_bound;
			HalfedgeLayer* hl_begin = (i == 0) ? cut.front()->hl : cut.back()->hl;
			HalfedgeLayer* hl_transfer;
			while (face_flag[hl_begin->left] || !face_flag[hl_begin->oppo->left])
			{
				hl_begin = hl_begin->prev->oppo;
			}
			hl_transfer = hl_begin;
			do
			{
				HalfedgeLayer* ht = hl_transfer->oppo;
				do
				{
					if (!face_flag[ht->left])
						break;
					ht = ht->prev->oppo;
				} while (true);
				hl_transfer = ht;
				one_bound.push_back(ht);
			} while (hl_transfer != hl_begin);
			bounds.push_back(std::move(one_bound));
		}
	}

	void cylinder::parametrize(M4 &m4, const Eigen::Matrix3Xd& normal)
	{
		//标记与cut相关的顶点和面，计算装配矩阵需要的数据
		int nvl = vertice_flag.size();
		int nfl = face_flag.size();
		BoolVector cutv_flag(nvl, false);
		BoolVector cutf_flag(nfl, false);
		std::vector<std::map<VertexLayer*, std::pair<bool, Eigen::Vector3d>>> info(nfl);
		{
			for (auto c : cut)
				cutv_flag[c->id] = true;
			HalfedgeLayer* hl = m4.find_halfedge_layer(cut[0], cut[1]);
			FaceLayer* fl = &m4.facelayers[hl->left];
			int flid = fl->id;
			int fid = fl->f.idx();
			Mesh* mesh = m4.mesh;
			auto calc_vector = [&](bool flag)
			{
				double inv_area = 1.0 / (2 * mesh->calc_face_area(fl->f));
				auto vi = hl->to;
				auto ev = mesh->calc_edge_vector(mesh->prev_halfedge_handle(hl->h));
				info[flid].insert(std::make_pair(&m4.verticelayers[vi],
					std::make_pair(flag && cutv_flag[vi], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
				vi = hl->next->to;
				ev = mesh->calc_edge_vector(hl->h);
				info[flid].insert(std::make_pair(&m4.verticelayers[vi],
					std::make_pair(flag && cutv_flag[vi], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
				vi = hl->from;
				ev = mesh->calc_edge_vector(mesh->next_halfedge_handle(hl->h));
				info[flid].insert(std::make_pair(&m4.verticelayers[vi],
					std::make_pair(flag && cutv_flag[vi], normal.col(fid).cross(Eigen::Vector3d(ev[0], ev[1], ev[2])) * inv_area)));
			};

			while (face_flag[flid])
			{
				cutf_flag[flid] = true;
				calc_vector(true);
				hl = hl->prev->oppo;
				fl = &m4.facelayers[hl->left];
				flid = fl->id;
				fid = fl->f.idx();
			}
			for (int i = 1; i < cut.size() - 1; ++i)
			{
				hl = m4.find_halfedge_layer(cut[i], cut[i + 1]);
				fl = &m4.facelayers[hl->left];
				flid = fl->id;
				fid = fl->f.idx();
				do
				{
					if (info[flid].empty() && face_flag[flid])
					{
						cutf_flag[flid] = true;
						calc_vector(true);
					}
					hl = hl->prev->oppo;
					fl = &m4.facelayers[hl->left];
					flid = fl->id;
					fid = fl->f.idx();
				} while (hl->to != cut[i - 1]->id);
			}
			hl = m4.find_halfedge_layer(cut[cut.size() - 2], cut.back())->next;
			fl = &m4.facelayers[hl->left];
			flid = fl->id;
			fid = fl->f.idx();
			while (face_flag[flid])
			{
				if (info[flid].empty())
				{
					cutf_flag[flid] = true;
					calc_vector(true);
				}
				hl = hl->oppo->next;
				fl = &m4.facelayers[hl->left];
				flid = fl->id;
				fid = fl->f.idx();
			}

			for (auto fa : faces)
			{
				if (cutf_flag[fa->id])
					continue;
				hl = fa->hl;
				fl = fa;
				flid = fl->id;
				fid = fl->f.idx();
				calc_vector(false);
			}
		}

		vidmap.resize(nvl, -1);
		int count = 0;
		for (auto vl : vertices)
			vidmap[vl->id] = count++;

		int vl_size = vertices.size();
		std::vector<Eigen::Triplet<double>> triple;
		std::vector<double> w(nvl, 0);
		uv[0].resize(vl_size); uv[0].setZero();
		uv[1].resize(vl_size); uv[1].setZero();

		auto &crossfield = m4.cf->crossfield;
		for (int i = 1; i < vl_size; ++i)
		{
			VertexLayer* vl = vertices[i];
			int vlid = vl->id;
			int vm = vidmap[vlid];
			HalfedgeLayer* hl_begin = vl->hl;
			HalfedgeLayer* hl_transfer = hl_begin;
			do
			{
				int vfl_id = hl_transfer->left;
				if (face_flag[vfl_id])
				{
					Eigen::Vector3d& R = info[vfl_id][vl].second;
					for (const auto& f_info : info[vfl_id])
					{
						double dot_ = R.dot(f_info.second.second);
						int fvlid = f_info.first->id;
						if (f_info.second.first)
							uv[0](vm) -= dot_;
						if (fvlid != vertices[0]->id)
							w[fvlid] += dot_;
					}
					uv[0](vm) += crossfield.col(vfl_id).dot(R);
					uv[1](vm) += crossfield.col(vfl_id + ((vfl_id % 4) == 3 ? -3 : 1)).dot(R);
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			triple.emplace_back(vm - 1, vm - 1, w[vlid]);
			w[vlid] = 0;

			hl_transfer = hl_begin;
			do
			{
				int vvid = hl_transfer->to;
				if (vertice_flag[vvid] && vvid != vertices[0]->id)
				{
					triple.emplace_back(vm - 1, vidmap[vvid] - 1, w[vvid]);
					w[vvid] = 0;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}
		static int es = 0;
		++es;
		if (es == 4)
		{
			int p = 0;
		}
		Eigen::SparseMatrix<double> A(vl_size - 1, vl_size - 1);
		A.setFromTriplets(triple.begin(), triple.end());
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		uv[0].tail(vl_size - 1) = solver.solve(uv[0].tail(vl_size - 1));
		uv[1].tail(vl_size - 1) = solver.solve(uv[1].tail(vl_size - 1));
	}
}