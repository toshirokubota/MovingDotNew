#pragma once
#include <CoreParticle.h>

namespace TK
{
	class MovingDot {
	public:
		MovingDot(CoreParticle* p)
		{
			this->x = p->x;
			this->y = p->y;
			this->z = p->z;
			this->t = p->t;
			dx = dy = dz = 0.0f;
			this->p = p;
			match = NULL;
		}
		int x;
		int y;
		int z;
		int t;
		float dx;
		float dy;
		float dz;
		CoreParticle* p;
		MovingDot* match;
	};
};
