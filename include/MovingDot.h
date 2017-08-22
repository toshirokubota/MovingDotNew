#pragma once
#include <CoreParticle.h>

namespace TK
{
	class MovingDot {
	public:
		static int _id;
		MovingDot(CoreParticle* p)
		{
			this->x = p->x;
			this->y = p->y;
			this->z = p->z;
			this->t = p->t;
			dx = dy = dz = 0.0f;
			this->p = p;
			label = 0;
			match = NULL;
			id = _id++;
		}
		int x;
		int y;
		int z;
		int t;
		int id;
		float dx;
		float dy;
		float dz;
		CoreParticle* p;
		MovingDot* match;
		int label;
	};
};
