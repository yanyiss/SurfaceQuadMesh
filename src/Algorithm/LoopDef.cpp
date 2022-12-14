#include "LoopDef.h"
namespace LoopGen
{
	void cylinder_set::push_back(cylinder& cy)
	{
		if (cylinders.empty())
		{
			cylinders.push_back(std::move(cy));
			return;
		}
		bool if_overlap = false;
		for (auto v : cy.vertices)
		{
			if (!vertex_cylinder_map[v.idx()].empty())
			{
				if_overlap = true;
				break;
			}
		}
		if (!if_overlap)
			return;

	}
}