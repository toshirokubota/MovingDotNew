#pragma once
#include <szParticleF.h>
#include <vector>
#include <set>
using namespace std;

class StationaryParticle;
class MovingParticle;
class PolygonFactory;

class Polygon
{
public:
	friend PolygonFactory;
	vector<CParticleF> project(float time);
	vector<CParticleF> original(); //returns locations of initial particles.
	bool operator==(const Polygon& poly) const
	{
		if (poly.particles.size() != particles.size()) return false;
		return sset == poly.sset;
	}

	int getId() const { return id; }
	float getCreatedTime() const { return creation_time; }
	int getNumParticles() const { return particles.size(); }
	vector<MovingParticle*> getParticles() const { return particles; }
	MovingParticle* getParticle(int i) { return (i < particles.size() && i >= 0 ? particles[i] : NULL); }
	bool contains(StationaryParticle* p) const { return p!=NULL && sset.find(p) != sset.end(); }
	bool contains(MovingParticle* p) const { return p != NULL && pset.find(p) != pset.end(); }
	set<MovingParticle*>& getParticleSet() { return pset; }
	set<StationaryParticle*>& getStationaryParticleSet() { return sset; }

	void print();

private:
	Polygon(vector<MovingParticle*>& vp, float creation_time, int t = 0);
	vector<MovingParticle*> particles;
	set<MovingParticle*> pset;
	set<StationaryParticle*> sset;
	float creation_time;
	int type;
	int id;
	static int _id;
};


class PolygonFactory
{
public:
	static PolygonFactory& getInstance()
	{
		static PolygonFactory* _instance = NULL;
		if (_instance == NULL)
		{
			_instance = new PolygonFactory();
		}
		return *_instance;
	}

	Polygon* makePolygon(vector<MovingParticle*>& vp, float ct)
	{
		Polygon* p = new Polygon(vp, ct);
		polygons.push_back(p);
		return p;
	}
	void clean()
	{
		for (int i = 0; i < polygons.size(); ++i)
		{
			delete polygons[i];
		}
		polygons.clear();
		Polygon::_id = 0;
	}

private:
	PolygonFactory()
	{
		Polygon::_id = 0;
	}
	~PolygonFactory()
	{
		clean();
	}
	PolygonFactory(PolygonFactory& f)
	{
	}
	PolygonFactory operator=(PolygonFactory& f)
	{
	}
	vector<Polygon*> polygons;
};
