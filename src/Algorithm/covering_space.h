#pragma once
#include "..\MeshViewer\MeshDefinition.h"
#include "..\Algorithm\crossField.h"
namespace LoopGen
{
	template <int layer>
	class CoveringSpaceDataStructure
	{
	public:
		CoveringSpaceDataStructure(Mesh* mesh_) : mesh(mesh_) { }
		~CoveringSpaceDataStructure() {}

	public:
		Mesh* mesh;
		crossField* cf;

		std::deque<bool> sing_flag;
		std::vector<std::vector<int>> vv;
		std::vector<std::vector<int>> fv;

		void set_cf(crossField* cf_) { cf = cf_; }
		void update()
		{
			int nv = mesh->n_vertices();
			int ne = mesh->n_edges();

			sing_flag.resize(nv, false);
			for (auto sing : cf->getSingularity())
				sing_flag[sing] = true;

			std::vector<int> h_shift(ne << 1);
			auto& matching = cf->getMatching();
			for (auto th = mesh->halfedges_begin(); th != mesh->halfedges_end(); ++th)
			{
				auto h = th.handle();
				auto fromv = mesh->from_vertex_handle(h);
				auto tov = mesh->to_vertex_handle(h);
				int fromid = fromv.idx();
				int toid = tov.idx();
				if ((sing_flag[fromid] || sing_flag[toid]) && fromid < toid)
					continue;

				int index = 0;
				auto h0 = mesh->voh_begin(fromv).handle();
				do
				{
					h0 = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h0));
					index += matching[h0.idx()]; index %= 4;
				} while (h0 != h);
				auto h1 = mesh->voh_begin(tov).handle();
				do
				{
					h1 = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(h1));
					index += 4 - matching[h1.idx()]; index %= 4;
				} while (h1 != h);

			}
		}
	};
	typedef CoveringSpaceDataStructure<2> M2;
	typedef CoveringSpaceDataStructure<4> M4;
}