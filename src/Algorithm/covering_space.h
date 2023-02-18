#pragma once
#include "..\MeshViewer\MeshDefinition.h"
#include "..\Algorithm\crossField.h"
#define YYSS_INFINITE 1.0e12
#define YYSS_FAIRLY_SMALL 1.0e-3
namespace LoopGen
{
	struct HalfedgeLayer 
	{
		HalfedgeHandle h; 
		//int layer;
		int id;
		HalfedgeLayer* prev; 
		HalfedgeLayer* next; 
		HalfedgeLayer* oppo; 
		int from; 
		int to; 
		int left; 
		HalfedgeLayer(){}
		HalfedgeLayer(HalfedgeHandle h_, /*int layer_, */int id_) : h(h_), /*layer(layer_), */id(id_), prev(nullptr), next(nullptr), oppo(nullptr) { }
		void set_info(HalfedgeLayer* prev_, HalfedgeLayer* next_, HalfedgeLayer* oppo_, int from_, int to_, int left_)
		{
			prev = prev_; next = next_; oppo = oppo_; from = from_; to = to_; left = left_;
		}
	};
	struct VertexLayer
	{ 
		VertexHandle v;
		//int layer;
		int id;
		HalfedgeLayer* hl; 
		VertexLayer() {} 
		VertexLayer(VertexHandle v_, /*int layer_, */int id_) : v(v_), /*layer(layer_), */id(id_), hl(nullptr) {}
	};
	struct FaceLayer
	{
		FaceHandle f; 
		//int layer;
		int id;
		HalfedgeLayer* hl; 
		FaceLayer() {}
		FaceLayer(FaceHandle f_, /*int layer_, */int id_) :f(f_), /*layer(layer_), */id(id_), hl(nullptr) {}
	};

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

		BoolVector sing_flag;
		std::vector<int> verticemap;
		std::vector<VertexLayer> verticelayers;
		std::vector<HalfedgeLayer> halfedgelayers;
		std::vector<FaceLayer> facelayers;
		Eigen::VectorXd weight;
		//std::vector<int> m4_to_m2id;
		//std::vector<int> m2_to_m4id;

		void set_base(Mesh* mesh_, crossField* cf_) { mesh = mesh_; cf = cf_; }
		int calc_shift(VertexHandle v, HalfedgeHandle h)
		{
			HalfedgeHandle h_transfer = mesh->voh_begin(v).handle();
			int index = 0;
			auto &matching = cf->getMatching();
			while (h_transfer != h)
			{
				h_transfer = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h_transfer));
				index += matching[h_transfer.idx()];
			}
			return index % 4;
		}
		void update()
		{
			int nv = mesh->n_vertices();
			int nh = mesh->n_halfedges();
			int nf = mesh->n_faces();

			sing_flag.resize(nv, false);
			for (auto sing : cf->getSingularity())
				sing_flag[sing] = true;

			verticemap.reserve(nv);
			{
				//确定顶点列表
				verticelayers.reserve(layer * nv - (layer - 1) * cf->getSingularity().size());
				//m4id.reserve(verticelayers.size());
				//m2id.reserve(layer / 2 * nv - (layer / 2 - 1)*cf->getSingularity().size());
				int idcount = 0;
				for (auto tv : mesh->vertices())
				{
					verticemap.push_back(idcount);
					if (sing_flag[tv.idx()])
					{
						verticelayers.emplace_back(tv, /*0, */idcount);
						++idcount;
					}
					else
					{
						for (int i = 0; i < layer; ++i)
						{
							verticelayers.emplace_back(tv, /*i, */idcount);
							++idcount;
						}
					}
				}
				//确定半边列表
				halfedgelayers.reserve(layer * nh);
				idcount = 0;
				for (auto th : mesh->halfedges())
				{
					for (int i = 0; i < layer; ++i)
					{
						halfedgelayers.emplace_back(th, /*i, */idcount);
						++idcount;
					}
				}
				//确定面列表
				facelayers.reserve(layer * nf);
				idcount = 0;
				for (auto tf : mesh->faces())
				{
					for (int i = 0; i < layer; ++i)
					{
						facelayers.emplace_back(tf, /*i, */idcount);
						//facelayers.emplace_back(tf, i, idcount);
						//facelayers.emplace_back(tf, i, idcount);
						++idcount;
					}
				}
			}

			for (auto tf : mesh->faces())
			{
				int fid = tf.idx();

				std::vector<HalfedgeHandle> hhs(3);
				std::vector<VertexHandle> vhs(3);
				std::vector<int> vshift(3);
				bool flag[3];
				hhs[0] = mesh->fh_begin(tf).handle();
				hhs[1] = mesh->next_halfedge_handle(hhs[0]);
				hhs[2] = mesh->next_halfedge_handle(hhs[1]);
				for (int i = 0; i < 3; ++i)
				{
					vhs[i] = mesh->from_vertex_handle(hhs[i]);
					vshift[i] = calc_shift(vhs[i], hhs[i]);
					flag[i] = sing_flag[vhs[i].idx()];
				}
				

				std::vector<HalfedgeLayer*> hlx(3);
				std::vector<int> vmx(3);
				for (int j = 0; j < layer; ++j)
				{
					for (int k = 0; k < 3; ++k)
					{
						hlx[k] = &halfedgelayers[hhs[k].idx()*layer + j];
						vmx[k] = verticemap[vhs[k].idx()] + (flag[k] ? 0 : (j + vshift[k]) % layer);
					}
					int left_face = fid * layer + j;
					hlx[0]->set_info(hlx[2], hlx[1], nullptr, vmx[0], vmx[1], left_face);
					hlx[1]->set_info(hlx[0], hlx[2], nullptr, vmx[1], vmx[2], left_face);
					hlx[2]->set_info(hlx[1], hlx[0], nullptr, vmx[2], vmx[0], left_face);
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
					halfedgelayers[h1.idx()*layer + (j + shift) % layer].oppo = &halfedgelayers[h0.idx()*layer + j];
					halfedgelayers[h0.idx()*layer + j].oppo = &halfedgelayers[h1.idx()*layer + (j + shift) % layer];
				}
			}

			//为顶点和面设定其遍历的第一条半边
			for (auto &tv : verticelayers)
			{
				HalfedgeHandle h = mesh->voh_begin(tv.v).handle();
				for (int j = 0; j < layer; ++j)
				{
					if (halfedgelayers[h.idx() * layer + j].from == tv.id)
					{
						tv.hl = &halfedgelayers[h.idx() * layer + j];
						break;
					}
				}
			}
			for (auto &tf : facelayers)
			{
				HalfedgeHandle h = mesh->fh_begin(tf.f).handle();
				for (int j = 0; j < layer; ++j)
				{
					if (halfedgelayers[h.idx() * layer + j].left == tf.id)
					{
						tf.hl = &halfedgelayers[h.idx() * layer + j];
						break;
					}
				}
			}
		}
		HalfedgeLayer* find_halfedge_layer(VertexLayer* vl0, VertexLayer* vl1)
		{
			auto hl_begin = vl0->hl;
			auto hl_transfer = hl_begin;
			do
			{
				if (hl_transfer->to == vl1->id)
				{
					return hl_transfer;
				}
				hl_transfer = hl_transfer->prev->oppo;
			} while (hl_transfer != hl_begin);
			return nullptr;
		}
		VertexLayer* conj_vl(VertexLayer* vl, int layer_plus)
		{
			int id = verticemap[vl->v.idx()];
			return &verticelayers[id + (vl->id - id + layer_plus) % 4];
		}
		void set_weight(double alpha = 900)
		{
			auto &crossfield = cf->getCrossField();
			auto &position = cf->getPosition();
			auto &normal = cf->getNormal();
			weight.resize(halfedgelayers.size());
			double quarterPI = PI * 0.25;
			//for (auto te : mesh->edges())
			for (auto th : mesh->halfedges())
			{
				HalfedgeHandle h0 = th;
				for (int j = 0; j < layer; ++j)
				{
					HalfedgeLayer& hl = halfedgelayers[h0.idx() * layer + j];
					auto &fv = crossfield.col(hl.left);
					auto &gv = crossfield.col(hl.oppo->left);
					auto ev = position.col(mesh->to_vertex_handle(h0).idx()) - position.col(mesh->from_vertex_handle(h0).idx());
					double arc0 = atan2(ev.cross(fv).dot(normal.col(hl.left / layer)), ev.dot(fv));
					double arc1 = atan2(ev.cross(gv).dot(normal.col(hl.oppo->left / layer)), ev.dot(gv));
					double arc = fabs(atan2(sin(arc0) + sin(arc1), cos(arc0) + cos(arc1)));
					if (arc < quarterPI)
						weight(hl.id) = ev.norm() * sqrt(alpha*sin(arc)*sin(arc) + 1);
					else
						weight(hl.id) = YYSS_INFINITE;
				}
			}
		}
	};
	//typedef CoveringSpaceDataStructure<2> M2;
	typedef CoveringSpaceDataStructure<4> M4;
}