#pragma once
#include "..\MeshViewer\MeshDefinition.h"
#include "..\Algorithm\crossField.h"
namespace LoopGen
{
	template <int layer>
	class CoveringSpaceDataStructure
	{
	public:
		CoveringSpaceDataStructure(Mesh* mesh_) : base_mesh(mesh_) { }
		~CoveringSpaceDataStructure() {}

	public:
		Mesh* base_mesh;
		crossField* cf;

		std::vector<std::vector<int>> vv;
		std::vector<std::vector<int>> fv;

		void set_cf(crossField* cf_) { cf = cf_; }
		void update()
		{

		}
	};
	typedef CoveringSpaceDataStructure<2> M2;
	typedef CoveringSpaceDataStructure<4> M4;
}