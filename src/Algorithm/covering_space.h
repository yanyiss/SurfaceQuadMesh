#pragma once
#include "..\MeshViewer\MeshDefinition.h"
#include "..\Algorithm\crossField.h"
namespace LoopGen
{
	struct VertexLayer { VertexHandle v; int id; VertexLayer() {} VertexLayer(VertexHandle v_, int id_) : v(v_), id(id_) {} };
	struct HalfedgeLayer 
	{
		HalfedgeHandle h; int id; HalfedgeLayer* prev; HalfedgeLayer* next; HalfedgeLayer* oppo; int from; int to; int left; 
		HalfedgeLayer(){}
		HalfedgeLayer(HalfedgeHandle h_, int id_) : h(h_), id(id_), prev(nullptr), next(nullptr), oppo(nullptr) { }
		void set_info(HalfedgeLayer* prev_, HalfedgeLayer* next_, HalfedgeLayer* oppo_, int from_, int to_, int left_)
		{
			prev = prev_; next = next_; oppo = oppo_; from = from_; to = to_; left = left_;
		}
	};
	struct FaceLayer { FaceHandle f; int id; FaceLayer() {} FaceLayer(FaceHandle f_, int id_) :f(f_), id(id_) {} };

	template <int layer>
	class CoveringSpaceDataStructure
	{
	public:
		CoveringSpaceDataStructure(){}
		~CoveringSpaceDataStructure() {}

	public:
		Mesh* mesh;
		crossField* cf;
		//int layer = 4;

		std::deque<bool> sing_flag;
		std::vector<VertexLayer> vertices;
		std::vector<HalfedgeLayer> halfedges;
		std::vector<FaceLayer> faces;

		void set_base(Mesh* mesh_, crossField* cf_) { mesh = mesh_; cf = cf_; }
		int calc_shift(VertexHandle v, HalfedgeHandle h)
		{
			HalfedgeHandle h_begin = mesh->voh_begin(v);
			int index = 0;
			auto &matching = cf->getMatching();
			while (h_begin != h)
			{
				h_begin = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h_begin));
				index += matching[h_begin.idx()];
			}
			return index % 4;
		}
		void update()
		{
			int nv = mesh->n_vertices();
			int ne = mesh->n_edges();
			int nf = mesh->n_faces();

			sing_flag.resize(nv, false);
			for (auto sing : cf->getSingularity())
				sing_flag[sing] = true;

			std::vector<int> verticemap; verticemap.reserve(nv);
			{
				//确定顶点列表
				vertices.reserve(layer * nv - (layer - 1) * cf->getSingularity().size());
				int idcount = 0;
				for (auto tv : mesh->vertices())
				{
					verticemap.push_back(idcount);
					if (sing_flag[tv.idx()])
					{
						vertices.emplace_back(tv, idcount);
						++idcount;
					}
					else
					{
						for (int i = 0; i < layer; ++i)
						{
							vertices.emplace_back(tv, idcount);
							++idcount;
						}
					}
				}
				//确定半边列表
				halfedges.reserve(2 * layer * ne);
				idcount = 0;
				for (auto th : mesh->halfedges())
				{
					for (int i = 0; i < layer; ++i)
					{
						halfedges.emplace_back(th, idcount);
						++idcount;
					}
				}
				//确定面列表
				faces.reserve(layer * nf);
				idcount = 0;
				for (auto tf : mesh->faces())
				{
					for (int i = 0; i < layer; ++i)
					{
						faces.emplace_back(tf, idcount);
						++idcount;
					}
				}
			}

			for (auto tf : mesh->faces())
			{
				int fid = tf.idx();
				std::vector<VertexHandle> vhs; vhs.reserve(3);
				std::deque<bool> flag;
				for (auto tfv : mesh->fv_range(tf))
				{
					vhs.push_back(tfv);
					if (sing_flag[tfv.idx()])
						flag.push_back(true);
					else
						flag.push_back(false);
				}
				std::swap(vhs[0], vhs[1]);
				std::swap(flag[0], flag[1]);
				
				std::vector<HalfedgeHandle> hhs; hhs.reserve(3);
				hhs.push_back(mesh->find_halfedge(vhs[0], vhs[1]));
				hhs.push_back(mesh->find_halfedge(vhs[1], vhs[2]));
				hhs.push_back(mesh->find_halfedge(vhs[2], vhs[0]));
				std::vector<int> vshift; hhs.reserve(3);
				vshift.push_back(calc_shift(vhs[0], hhs[0]));
				vshift.push_back(calc_shift(vhs[1], hhs[1]));
				vshift.push_back(calc_shift(vhs[2], hhs[2]));

				for (int j = 0; j < layer; ++j)
				{
					halfedges[hhs[0].idx() * layer + j].set_info(&halfedges[hhs[2].idx() * layer + j],
						&halfedges[hhs[1].idx() * layer + j], nullptr, verticemap[vhs[0].idx()] + flag[0] ? (j + vshift[0]) % layer : 0,
						verticemap[vhs[1].idx()] + flag[1] ? (j + vshift[1]) % layer : 0, fid * layer + j);
					halfedges[hhs[1].idx() * layer + j].set_info(&halfedges[hhs[0].idx() * layer + j],
						&halfedges[hhs[2].idx() * layer + j], nullptr, verticemap[vhs[1].idx()] + flag[1] ? (j + vshift[1]) % layer : 0,
						verticemap[vhs[2].idx()] + flag[2] ? (j + vshift[2]) % layer : 0, fid * layer + j);
					halfedges[hhs[2].idx() * layer + j].set_info(&halfedges[hhs[1].idx() * layer + j],
						&halfedges[hhs[0].idx() * layer + j], nullptr, verticemap[vhs[2].idx()] + flag[2] ? (j + vshift[2]) % layer : 0,
						verticemap[vhs[0].idx()] + flag[0] ? (j + vshift[0]) % layer : 0, fid * layer + j);
				}
			}

			//确定相反半边的关系
			auto &matching = cf->getMatching();
			for (auto te : mesh->edges())
			{
				VertexHandle v0 = te.v0();
				VertexHandle v1 = te.v1();
				HalfedgeHandle h0 = mesh->find_halfedge(v0, v1);
				HalfedgeHandle h1 = mesh->find_halfedge(v1, v0);
				int shift = matching[h0.idx()];
				for (int j = 0; j < layer; ++j)
				{
					halfedges[h1.idx()*layer + (j + shift) % layer].oppo = &halfedges[h0.idx()*layer + j];
					halfedges[h0.idx()*layer + j].oppo = &halfedges[h1.idx()*layer + (j + shift) % layer];
				}
			}

		}
	};
	typedef CoveringSpaceDataStructure<2> M2;
	typedef CoveringSpaceDataStructure<4> M4;
}