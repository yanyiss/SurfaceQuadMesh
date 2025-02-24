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

		BoolVector sing_flag;
		std::vector<VertexLayer> verticelayers;
		std::vector<HalfedgeLayer> halfedgelayers;
		std::vector<FaceLayer> facelayers;
		Eigen::VectorXd weight;
		Eigen::VectorXd weight_2ring;

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
		void init()
		{
			//确定顶点列表
			verticelayers.clear();
			verticelayers.reserve(layer*mesh->n_vertices());
			int idcount = 0;
			for (auto tv : mesh->vertices())
			{
				for (int i = 0; i < layer; ++i)
				{
					verticelayers.emplace_back(tv, idcount);
					++idcount;
				}
			}
			//确定半边列表
			halfedgelayers.clear();
			halfedgelayers.reserve(layer * mesh->n_halfedges());
			idcount = 0;
			for (auto th : mesh->halfedges())
			{
				for (int i = 0; i < layer; ++i)
				{
					halfedgelayers.emplace_back(th, idcount);
					++idcount;
				}
			}
			//确定面列表
			facelayers.clear();
			facelayers.reserve(layer * mesh->n_faces());
			idcount = 0;
			for (auto tf : mesh->faces())
			{
				for (int i = 0; i < layer; ++i)
				{
					facelayers.emplace_back(tf, idcount);
					++idcount;
				}
			}
		}
		void update()
		{
			sing_flag.resize(mesh->n_vertices(), false);
			for (auto &sing : cf->getSingularity())
				for (auto &s : sing)
					sing_flag[s] = true;

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
						vmx[k] = vhs[k].idx() * 4 + (flag[k] ? 0 : ((j + vshift[k]) % layer));
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
		VertexLayer* another_layer(VertexLayer* vl, int layer_plus)
		{
			//int id = verticemap[vl->v.idx()];
			int id = vl->v.idx() * 4;
			return &verticelayers[id + (vl->id - id + layer_plus) % 4];
		}
		FaceLayer* another_layer(FaceLayer* fl, int layer_plus)
		{
			int id = fl->f.idx() * 4;
			return &facelayers[id + (fl->id - id + layer_plus) % 4];
		}
		HalfedgeLayer* another_layer(HalfedgeLayer* hl, int layer_plus)
		{
			int id = hl->h.idx() * 4;
			return &halfedgelayers[id + (hl->id - id + layer_plus) % 4];
		}
		void set_weight(double alpha = 900)
		{
			auto &crossfield = cf->getCrossField();
			auto &position = cf->getPosition();
			auto &normal = cf->getNormal();
			weight.resize(halfedgelayers.size());
			weight_2ring.resize(halfedgelayers.size());
			double quarterPI = PI * 0.25;
			double halfPI = PI * 0.5;
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
					
#if 0
					double arc = fabs(atan2(sin(arc0) + sin(arc1), cos(arc0) + cos(arc1)));
					if (arc < quarterPI)
						weight(hl.id) = ev.norm() * sqrt(alpha * sin(arc) * sin(arc) + 1);
					/*{
						weight(hl.id) = sin(arc);
						weight(hl.id) *= weight(hl.id);
					}*/
					/*else if (arc < halfPI)
						weight(hl.id) = YYSS_INFINITE / 2;*/
					else
						weight(hl.id) = YYSS_INFINITE;
#else
					//if (arc0 < 0) arc0 += 2 * PI;
					//if (arc1 < 0) arc1 += 2 * PI;
					double arc = (fabs(arc0) + fabs(arc1)) / 2;
					if (arc < quarterPI)
						weight(hl.id) = ev.norm() * sqrt(alpha * sin(arc) * sin(arc) + 1);
					//else if (arc < 1.57)
					//	weight(hl.id) = ev.norm() * sqrt(alpha * sin(arc) * sin(arc) + 1);
					else
						weight(hl.id) = YYSS_INFINITE;
#endif
				}
			}

			for (auto th : mesh->halfedges())
			{
				HalfedgeHandle h0 = th;
				for (int j = 0; j < layer; ++j)
				{
					HalfedgeLayer& hl = halfedgelayers[h0.idx() * layer + j];
					auto& fv = crossfield.col(hl.left);
					auto& gv = crossfield.col(hl.next->oppo->left);
					auto n0 = normal.col(hl.left / 4);
					auto n1 = normal.col(hl.next->oppo->left / 4);

					auto p0 = position.col(hl.from / 4);
					auto p1 = position.col(hl.to / 4);
					auto p2 = position.col(hl.next->to / 4);
					auto p3 = position.col(hl.next->oppo->next->to / 4);





					auto ev = position.col(mesh->to_vertex_handle(h0).idx()) - position.col(mesh->from_vertex_handle(h0).idx());
					double arc0 = atan2(ev.cross(fv).dot(normal.col(hl.left / layer)), ev.dot(fv));
					double arc1 = atan2(ev.cross(gv).dot(normal.col(hl.oppo->left / layer)), ev.dot(gv));

					double arc = (fabs(arc0) + fabs(arc1)) / 2;
					if (arc < quarterPI)
						weight(hl.id) = ev.norm() * sqrt(alpha * sin(arc) * sin(arc) + 1);
					else
						weight(hl.id) = YYSS_INFINITE;
				}
			}
		}
	};
	//typedef CoveringSpaceDataStructure<2> M2;
	typedef CoveringSpaceDataStructure<4> M4;

#define vhMarkCirculator(hl_begin, hl_transfer, codeBody)\
do\
{\
##codeBody##\
##hl_transfer## = ##hl_transfer##->prev->oppo;\
}while(##hl_transfer## != ##hl_begin##);

#define fhMarkCirculator(hl_begin, hl_transfer, codeBody)\
do\
{\
##codeBody##\
##hl_transfer## = ##hl_transfer##->next;\
}while(##hl_transfer## != ##hl_begin##);

#define vhCirculator(vl, codeBody)\
{\
HalfedgeLayer* hl_begin = ##vl##->hl;\
HalfedgeLayer* hl_transfer = hl_begin;\
do\
{\
##codeBody##\
hl_transfer = hl_transfer->prev->oppo;\
}while(hl_transfer != hl_begin);\
}

#define fhCirculator(fl, codeBody)\
{\
HalfedgeLayer* hl_begin = ##fl##->hl;\
HalfedgeLayer* hl_transfer = hl_begin;\
do\
{\
##codeBody##\
hl_transfer = hl_transfer->next;\
}while(hl_transfer != hl_begin);\
}

}