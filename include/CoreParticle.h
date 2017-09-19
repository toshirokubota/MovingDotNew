#pragma once
#include <vector>
#include <set>
using namespace std;
#include <Graph.h>

struct CoreParticle
{
	CoreParticle(int x0 = 0, int y0 = 0, int z0 = 0, int t0 = 0, int gen0 = 0, int lb = 0, int id0 = 0)
	{
		x = x0;
		y = y0;
		gen = gen0;
		label = lb;
		id = id0;
		z = z0;
		t = t0;
		value = 0;
		dval = 0;
		src = NULL;
		//selected = false;
		core = -1; //this will be set either 0, 1, or 2.
		vertex = NULL;
		//saliency = std::numeric_limits<float>::quiet_NaN();
	}
	int x;
	int y;
	int z;
	int t;
	int label;
	int gen; //generation
	int id;
	float value; //generic value
	float dval; //distance value
	//float saliency; //new addition 08/09/2017
	int core;

	set<CoreParticle*> sources;
	set<CoreParticle*> ascendents;
	set<CoreParticle*> descendents;
	vector<CoreParticle*> neighbors;
	vector<CoreParticle*> neighbors8;
	//CoreParticle* pi; //path from the core
	CoreParticle* src; //strong core of this particle
	TK::Vertex<CoreParticle*>* vertex;
};

struct CoreParticleFactory
{
public:
	static CoreParticleFactory& getInstance()
	{
		static CoreParticleFactory instance;
		return instance;
	}
	CoreParticle* makeParticle(int x = 0, int y = 0, int z = 0, int t = 0, int gen = 0, int lb = 0)
	{
		CoreParticle* particle = new CoreParticle(x, y, z, t, gen, lb, _id++);
		particles.push_back(particle);
		return particle;
	}
	int deleteParticle(set<CoreParticle*>& toremove)
	{
		int count = 0;
		for (int i = particles.size()-1; i>=0; --i)
		{
			if (toremove.find(particles[i]) != toremove.end())
			{
				delete particles[i];
				particles.erase(particles.begin() + i);
				count++;
			}
		}
		return count;
	}
	void clean()
	{
		for (int i = 0; i < particles.size(); ++i)
		{
			delete particles[i];
		}
		particles.clear();
		_id = 0;
	}
	vector<CoreParticle*> particles;
private:
	int _id;
	CoreParticleFactory()
	{
		_id = 0;
	}
	~CoreParticleFactory()
	{
		clean();
		_id = 0;
	}
	CoreParticleFactory(CoreParticleFactory& f) {}
	CoreParticleFactory operator=(CoreParticleFactory& f) {}
};
