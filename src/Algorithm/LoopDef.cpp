#include "LoopDef.h"
namespace LoopGen
{
	void region::set_face(M4 &m4)
	{
		face_flag.resize(m4.facelayers.size(), false);
		for (auto vl : vertices)
		{
			HalfedgeLayer* hl_begin = vl->hl;
			HalfedgeLayer* hl_transfer = hl_begin;
			do
			{
				if (!face_flag[hl_transfer->left])
				{
					if (vertice_flag[hl_transfer->to] && vertice_flag[hl_transfer->next->to])
					{
						faces.push_back(&m4.facelayers[hl_transfer->left]);
						face_flag[hl_transfer->left] = true;
					}
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
		}
	}

	cylinder::cylinder(cylinder &&cy)
	{
		id = cy.id; cy.id = -1;
		vertices = std::move(cy.vertices);
		faces = std::move(cy.faces);
		cut = std::move(cy.cut);
		bounds = std::move(cy.bounds);
		vertice_flag = std::move(cy.vertice_flag);
		face_flag = std::move(cy.face_flag);
		vidmap = std::move(cy.vidmap);
		uv[0] = std::move(cy.uv[0]); uv[1] = std::move(cy.uv[1]);
		handle_to_layer = std::move(cy.handle_to_layer);
	}

	cylinder& cylinder::operator=(cylinder&& cy)
	{
		id = cy.id; cy.id = -1;
		vertices = std::move(cy.vertices);
		faces = std::move(cy.faces);
		cut = std::move(cy.cut);
		bounds = std::move(cy.bounds);
		vertice_flag = std::move(cy.vertice_flag);
		face_flag = std::move(cy.face_flag);
		vidmap = std::move(cy.vidmap);
		uv[0] = std::move(cy.uv[0]); uv[1] = std::move(cy.uv[1]);
		handle_to_layer = std::move(cy.handle_to_layer);
		return *this;
	}
	
	void cylinder::set_bound()
	{
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

		Eigen::SparseMatrix<double> A(vl_size - 1, vl_size - 1);
		A.setFromTriplets(triple.begin(), triple.end());
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		uv[0].tail(vl_size - 1) = solver.solve(uv[0].tail(vl_size - 1));
		uv[1].tail(vl_size - 1) = solver.solve(uv[1].tail(vl_size - 1));
	}

	void cylinder_set::ProcessOverlap(M4 &m4)
	{
		if (cylinders.size() < 2)
			return;
		int nvl = m4.verticelayers.size();
		std::vector<std::vector<int>> region_index(nvl);
		for (auto &cy : cylinders)
		{
			for (auto vl : cy.vertices)
			{
				region_index[vl->id].push_back(cy.id);
				region_index[m4.another_layer(vl, 2)->id].push_back(cy.id);
			}
		}
		for (int i = 0; i < nvl; ++i)
		{
			int ris = region_index[i].size();
			while (ris > 1)
			{
				int parent = region_index[i][ris - 2];
				int child = region_index[i][ris - 1];
				auto &pc = cylinders[parent];
				auto &cc = cylinders[child];
				if (pc.vertice_flag[i] && cc.vertice_flag[i])
				{
					for (auto vl : cc.vertices)
					{
						if (pc.vertice_flag[vl->id])
							continue;
						pc.vertices.push_back(vl);
						pc.vertice_flag[vl->id] = true;
					}
				}
				else
				{
					VertexLayer* oppo_vl = nullptr;
					for (auto vl : cc.vertices)
					{
						oppo_vl = m4.another_layer(vl, 2);
						if (pc.vertice_flag[oppo_vl->id])
							continue;
						pc.vertices.push_back(oppo_vl);
						pc.vertice_flag[oppo_vl->id] = true;
					}
				}
				cc.id = -1;
				for (int j = 0; j < nvl; ++j)
				{
					if (region_index[j].empty())
						continue;
					if (region_index[j].back() == child)
						region_index[j].pop_back();
				}
				ris = region_index[i].size();
			}
		}

		std::vector<cylinder> new_cylinders;
		for (auto &cy : cylinders)
		{
			if (cy.id == -1)
				continue;
			new_cylinders.push_back(std::move(cy));
			new_cylinders.back().id = new_cylinders.size() - 1;
		}
		cylinders = std::move(new_cylinders);
	}

	disk::disk(disk &&dk)
	{
		id = dk.id; dk.id = -1;
		vertices = std::move(dk.vertices);
		faces = std::move(dk.faces);
		vertice_flag = std::move(dk.vertice_flag);
		face_flag = std::move(dk.face_flag);
		bounds = std::move(dk.bounds);
	}

	disk& disk::operator=(disk&& dk)
	{
		id = dk.id; dk.id = -1;
		vertices = std::move(dk.vertices);
		faces = std::move(dk.faces);
		vertice_flag = std::move(dk.vertice_flag);
		face_flag = std::move(dk.face_flag);
		bounds = std::move(dk.bounds);
		return *this;
	}

	void disk::set_bound()
	{
		bounds.resize(2);
		
	}

	void disk_set::ProcessOverlap(M4 &m4)
	{
		if (disks.size() < 2)
			return;
		int nvl = m4.verticelayers.size();
		std::vector<std::vector<int>> region_index(nvl);
		for (auto &dk : disks)
		{
			for (auto vl : dk.vertices)
			{
				region_index[vl->id].push_back(dk.id);
				region_index[m4.another_layer(vl, 2)->id].push_back(dk.id);
			}
		}
		for (int i = 0; i < nvl; ++i)
		{
			int ris = region_index[i].size();
			while (ris > 1)
			{
				int parent = region_index[i][ris - 2];
				int child = region_index[i][ris - 1];
				auto &pd = disks[parent];
				auto &cd = disks[child];
				if (pd.faces.size() < cd.faces.size())
				{
					disk dk = std::move(pd);
					pd = std::move(cd);
					cd = std::move(dk);
				}
				cd.id = -1;
				for (int j = 0; j < nvl; ++j)
				{
					if (region_index[j].empty())
						continue;
					if (region_index[j].back() == child)
						region_index[j].pop_back();
				}
				ris = region_index[i].size();
			}
		}
		std::vector<disk> new_disks;
		for (auto &dk : disks)
		{
			if (dk.id == -1)
				continue;
			new_disks.push_back(std::move(dk));
			new_disks.back().id = new_disks.size() - 1;
		}
		disks = std::move(new_disks);
	}
}