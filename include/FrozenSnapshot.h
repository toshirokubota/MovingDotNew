#pragma once
#include <vector>
using namespace std;
#include <szParticleF.h>
#include <mex.h>

class FrozenSnapshot
{
public:
	FrozenSnapshot(vector<CParticleF>& particles)
	{
		this->particles = particles;
	}
	static mxArray* StoreSnapshot(FrozenSnapshot& snapshot);
	static mxArray* StoreSnapshots(vector<FrozenSnapshot>& snapshot);

	vector<CParticleF> getParticles() const { return particles; }
protected:
	vector<CParticleF> particles;
};
