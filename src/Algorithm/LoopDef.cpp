#include "LoopDef.h"
namespace LoopGen
{
	void cylinder::set_bound()
	{
		for (int i = 0; i < 2; ++i)
		{
			std::vector<HalfedgeLayer*> one_bound;
			HalfedgeLayer* hl_begin = (i == 0) ? cut.front()->hl : cut.back()->hl;
			HalfedgeLayer* hl_transfer;
			while (!face_flag[hl_begin->left])
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
}