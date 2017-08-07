#pragma once
#include <ParticleSimulator.h>

struct FitnessStruct
{
	FitnessStruct()
	{
		value = -std::numeric_limits<float>::infinity();
		coverage = 0;
		fitness = 0;
		bleft = true;
		f = NULL;
		p = NULL;
	}
	MovingParticle* p;
	vector<MovingParticle*> area;
	float fitness;
	float coverage;
	float value;
	bool bleft;
	MovingParticle* f;
};

class ParticleSimulatorGreedy : public ParticleSimulator
{
public:
	ParticleSimulatorGreedy() : ParticleSimulator()
	{
	}
	virtual bool Simulate(float endtime = 10.0f, float delta = 0.1f, bool bdebug = false);
	float _fitnessMeasure(vector<CParticleF>& vp, float scale);
	FitnessStruct findNextEventGreedy();
	pair<MovingParticle*, MovingParticle*> applyEventGreedy(EventStruct ev);
	virtual FitnessStruct computeFitness(EventStruct ev, bool left);
	MovingParticle* pickCollidingParticle(EventStruct ev);
	float _getCoverage(vector<MovingParticle*>& area);
	vector<Polygon*> chosen;
	set<StationaryParticle*> covered;
	static CParticleF makeSplitParticle(EventStruct ev);
	static vector<MovingParticle*> extractSimplePath(MovingParticle* p0, MovingParticle* pend);
};