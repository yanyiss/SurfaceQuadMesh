#include "LoopDef.h"
namespace LoopGen
{
	region_set rset = region_set();

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

	disk::disk(disk &&dk)
	{
		id = dk.id; dk.id = -1;
		vertices = std::move(dk.vertices);
		faces = std::move(dk.faces);
		vertice_flag = std::move(dk.vertice_flag);
		face_flag = std::move(dk.face_flag);
		bounds = std::move(dk.bounds);
		handle_to_layer = std::move(dk.handle_to_layer);
	}

	disk& disk::operator=(disk&& dk)
	{
		id = dk.id; dk.id = -1;
		vertices = std::move(dk.vertices);
		faces = std::move(dk.faces);
		vertice_flag = std::move(dk.vertice_flag);
		face_flag = std::move(dk.face_flag);
		bounds = std::move(dk.bounds);
		handle_to_layer = std::move(dk.handle_to_layer);
		return *this;
	}

	void disk::set_bound(M4 &m4)
	{
		bounds.resize(2);
		HalfedgeLayer* hl;
		int shift;

		shift = to_bound.second == 0 ? 1 : 3;
		for (auto hl_transfer : rset.regions[to_bound.first]->bounds[to_bound.second])
		{
			if (vertice_flag[m4.another_layer(hl_transfer, shift)->to]
				&& !vertice_flag[m4.another_layer(hl_transfer, shift)->from])
			{
				hl = m4.another_layer(hl_transfer, shift)->oppo;
				break;
			}
		}
		while (true)
		{
			if (!face_flag[hl->left] && face_flag[hl->oppo->left])
				break;
			hl = hl->oppo->next;
		}
		shift = from_bound.second == 0 ? 1 : 3;
		auto &fvf = rset.regions[from_bound.first]->vertice_flag;
		while (!fvf[m4.another_layer(hl, shift)->to])
		{
			bounds[0].push_back(hl);
			hl = hl->oppo;
			while (true)
			{
				if (!face_flag[hl->left] && face_flag[hl->oppo->left])
					break;
				hl = hl->oppo->next;
			}
		}
		bounds[0].push_back(hl);


		shift = from_bound.second == 0 ? 3 : 1;
		for (auto hl_transfer : rset.regions[from_bound.first]->bounds[from_bound.second])
		{
			if (vertice_flag[m4.another_layer(hl_transfer, shift)->to]
				&& !vertice_flag[m4.another_layer(hl_transfer, shift)->from])
			{
				hl = m4.another_layer(hl_transfer, shift)->oppo;
				break;
			}
		}
		while (true)
		{
			if (!face_flag[hl->left] && face_flag[hl->oppo->left])
				break;
			hl = hl->oppo->next;
		}
		shift = to_bound.second == 0 ? 3 : 1;
		auto &tvf = rset.regions[to_bound.first]->vertice_flag;
		while (!tvf[m4.another_layer(hl, shift)->to])
		{
			bounds[1].push_back(hl);
			hl = hl->oppo;
			while (true)
			{
				if (!face_flag[hl->left] && face_flag[hl->oppo->left])
					break;
				hl = hl->oppo->next;
			}
		}
		bounds[1].push_back(hl);
	}

	void region_set::ProcessOverlap(M4 &m4, int type)
	{
		if (type == 0)
		{
			if (regions.size() < 2)
				return;
			int nvl = m4.verticelayers.size();
			std::vector<std::vector<int>> region_index(nvl);
			for (auto &rg : regions)
			{
				for (auto vl : rg->vertices)
				{
					region_index[vl->id].push_back(rg->id);
					region_index[m4.another_layer(vl, 2)->id].push_back(rg->id);
				}
			}
			for (int i = 0; i < nvl; ++i)
			{
				int ris = region_index[i].size();
				while (ris > 1)
				{
					int parent = region_index[i][ris - 2];
					int child = region_index[i][ris - 1];
					region* pc = regions[parent];
					region* cc = regions[child];
					if (pc->vertice_flag[i] && cc->vertice_flag[i])
					{
						for (auto vl : cc->vertices)
						{
							if (pc->vertice_flag[vl->id])
								continue;
							pc->vertices.push_back(vl);
							pc->vertice_flag[vl->id] = true;
						}
					}
					else
					{
						VertexLayer* oppo_vl = nullptr;
						for (auto vl : cc->vertices)
						{
							oppo_vl = m4.another_layer(vl, 2);
							if (pc->vertice_flag[oppo_vl->id])
								continue;
							pc->vertices.push_back(oppo_vl);
							pc->vertice_flag[oppo_vl->id] = true;
						}
					}
					cc->id = -1;
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
			for (auto rg : regions)
			{
				if (rg->id == -1)
				{
					delete rg;
					rg = nullptr;
				}
			}
			regions.erase(std::remove_if(regions.begin(), regions.end(), [&](region* rhs) {return rhs->id == -1; }), regions.end());
			int count = 0;
			for (auto rg : regions)
				rg->id = count++;
		}
		else
		{
			if (regions.size() < disk_mark + 2)
				return;
			int nvl = m4.verticelayers.size();
			std::vector<std::vector<int>> region_index(nvl);
			for (auto &rg : regions)
			{
				if (rg->id < disk_mark)
					continue;
				for (auto vl : rg->vertices)
				{
					region_index[vl->id].push_back(rg->id);
					region_index[m4.another_layer(vl, 2)->id].push_back(rg->id);
				}
			}
			for (int i = 0; i < nvl; ++i)
			{
				int ris = region_index[i].size();
				while (ris > 1)
				{
					int parent = region_index[i][ris - 2];
					int child = region_index[i][ris - 1];
					region* pd = regions[parent];
					region* cd = regions[child];
					if (pd->faces.size() < cd->faces.size())
					{
						std::swap(regions[parent], regions[child]);
						std::swap(pd, cd);
					}
					cd->id = -1;
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
			for (auto rg : regions)
			{
				if (rg->id == -1)
				{
					delete rg;
					rg = nullptr;
				}
			}
			regions.erase(std::remove_if(regions.begin(), regions.end(), [&](region* rhs) {return rhs->id == -1; }), regions.end());
			int count = 0;
			for (auto rg : regions)
				rg->id = count++;
		}
	}

	SegmentTree::SegmentTree(PlaneLoop &pl, M4 &m4)
	{
		pos.reserve(pl.size());
		for (auto &poh : pl)
			pos.push_back(poh.point(m4));
	}

	double SegmentTree::closest_distance(OpenMesh::Vec3d &in)
	{
		OpenMesh::Vec3d p01, pin;
		double proj, norm;
		double dis = YYSS_INFINITE;
		for (int i = 0; i < pos.size() - 1; ++i)
		{
			auto &p0 = pos[i];
			auto &p1 = pos[i + 1];
			p01 = (p1 - p0).normalized();
			pin = in - p0;
			proj = pin.dot(p01);
			if (proj <= 0)
			{
				norm = pin.norm();
				if (norm < dis)
				{
					dis = norm;
				}
			}
			else
			{
				pin = in - p1;
				proj = pin.dot(-p01);
				if (proj <= 0)
				{
					norm = pin.norm();
					if (norm < dis)
					{
						dis = norm;
					}
				}
				else
				{
					pin = p1 + proj * (-p01);
					norm = (pin - in).norm();
					if (norm < dis)
					{
						dis = norm;
					}
				}
			}
		}
		return dis;
	}
}