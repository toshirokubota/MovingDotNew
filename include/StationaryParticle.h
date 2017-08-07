#pragma once
#include <vector>
#include <set>
#include <map>
using namespace std;
#include <szParticleF.h>

class MovingParticle;
class StationaryParticleFactory;

class StationaryParticle
{
public:
	friend StationaryParticleFactory; //so that the factory can access _id.

	void print(char* tab = "") const
	{
		printf("%s%d\t%f\t%f\n", tab, id, p.m_X, p.m_Y);
	}

	CParticleF getP() const { return p; }
	float getX() const { return p.m_X; }
	float getY() const { return p.m_Y; }
	int getId() const { return id; }
	void* link; //generic pointer used in ParticleSimulatorGreedy::extractClosedPolygon().
	int info; //generic data used in ParticleSimulatorGreedy::extractClosedPolygon().
	set<StationaryParticle*> neighbors;
	set<MovingParticle*> derived;
private:
	StationaryParticle(CParticleF& p = CParticleF(0, 0, 0), int i=-1)
	{
		this->p = p;
		id = i;
		link = NULL;
		info = -1;
	}

	CParticleF p;
	int id;
};

struct StationaryParticleFactory
{
public:
	static StationaryParticleFactory& getInstance()
	{
		static StationaryParticleFactory instance;
		return instance;
	}
	StationaryParticle* makeParticle(CParticleF& p)
	{
		StationaryParticle* particle = new StationaryParticle(p, _id++);
		particles.push_back(particle);
		return particle;
	}
	void clean()
	{
		for (int i = 0; i<particles.size(); ++i)
		{
			delete particles[i];
		}
		particles.clear();
		_id = 0;
	}
	vector<StationaryParticle*> particles;
private:
	int _id;
	StationaryParticleFactory()
	{
		_id = 0;
	}
	~StationaryParticleFactory()
	{
		clean();
	}
	StationaryParticleFactory(StationaryParticleFactory& f){}
	StationaryParticleFactory operator=(StationaryParticleFactory& f){}
};

