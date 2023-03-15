#include "AngleDefectMinimizer.h"
#include "../src/Toolbox/dprint.h"
namespace LoopGen
{
	void AngleDefectMinimizer::run()
	{
		auto ad = [&](VertexHandle v)
		{
			double angle_defect = 2 * M_PI;
			for (auto vih : mesh->vih_range(v))
				angle_defect -= mesh->calc_sector_angle(vih);
			return angle_defect;
		};
		for (auto tv = mesh->vertices_begin(); tv != mesh->vertices_end(); ++tv)
		{
			double angle_defect = ad(tv.handle());
			if (angle_defect > 0.2)
			{
				for (auto tvoh : mesh->voh_range(tv.handle()))
				{
					if (!mesh->is_flip_ok(tvoh.next().edge()))
						continue;
					if (ad(tvoh.to()) < -0.2 && ad(tvoh.next().to()) < -0.2 && ad(tvoh.next().opp().next().to()) > 0.1)
					{
						mesh->flip(tvoh.next().edge());
						break;
					}
				}
			}
		}
		mesh->garbage_collection();
	}
}