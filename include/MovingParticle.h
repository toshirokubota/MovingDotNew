#pragma once
#include <vector>
#include <set>
#include <map>
using namespace std;

#include <szParticleF.h>
#include <EventStruct.h>
#include <MiscGeometry.h>
#include <szMiscOperations.h>
#include <Graph.h>
#include <MovingFront.h>
#include <StationaryParticle.h>

enum MovingParticleType { Unknown, Initial, Regular, Merge, Split, Collide, Axis, Dummy };
struct ParticleFactory;
class ParticleSimulator;
class OffsetPolygonDAGBuilder;
class Polygon;

MovingParticleType int2ParticleType(int i);
const int ParticleDumpSize = 27;

class MovingParticle
{
public:
	static const int NumParents = 4;
	friend ParticleFactory; //so that the factory can access _id.
	friend ParticleSimulator;
	friend OffsetPolygonDAGBuilder;

	void update(float delta)
	{
		if (bActive && !bUnstable)
		{
			p = move(delta);
			time += delta;
		}
	}
	void print(char* tab = "") const;
	void printParentTree(char* tab = "") const;

	StationaryParticle* getInitParticle() const { return init_particle; }
	CParticleF getP0() const { return p0; }
	CParticleF getP() const { return p; }
	MovingFront getFront() const { return front; }
	MovingFront getRear() const { return rear; }
	void getVelocity(float& vx, float& vy) const { vx = v[0]; vy = v[1]; }
	MovingParticle* getNext() const { return next; }
	MovingParticle* getPrev() const { return prev; }
	int getId() const { return id; }
	EventStruct getEvent() const { return event; }
	float getTime() const { return time; }
	float getCreatedTime() const { return created; }
	MovingParticleType getType() const { return type; }
	MovingParticle* getParent(int j) const
	{
		if (j>=0 && j<NumParents) return this->parents[j];
		else return NULL;
	}
	Polygon* getPolygon() const { return polygon; }
	void setPolygon(Polygon* poly) { polygon = poly; }
	bool isActive() const { return bActive; }
	bool isReflex(); //check if it is concave (reflex) that allows splitting of a side
	bool isLeaf() const { return bLeaf; }
	bool isUnstable() const { return bUnstable; }
	void setEvent(EventStruct ev) {event = ev; }
	CParticleF move(float delta) const
	{
		return CParticleF(p.m_X + delta*v[0], p.m_Y + delta*v[1]);
	}
	CParticleF project(float time) const
	{
		float delta = time - created;
		return CParticleF(p0.m_X + delta*v[0], p0.m_Y + delta*v[1]);
	}
	bool updateEvent();
	bool initializeVelocity();
	bool calculateVelocityR();
	void clearVelocity()
	{
		v[0] = v[1] = 0.0f;
	}
	pair<MovingParticle*, MovingParticle*> applyEvent();
	vector<float> dump2vector(); //store the current state as a vector of float 
	EventStruct findNextEdgeEvent() const;
	EventStruct findNextSplitEvent() const;
	EventStruct findNextSplitEvent(const vector<MovingParticle*>& vp) const;

	static float frontPropAngle(MovingParticle* p, MovingParticle* q); //find the angle of propagating front p-q.
	static void setNeighbors(MovingParticle* p, MovingParticle* prev, MovingParticle* next)
	{
		setNeighbors(p, prev, next, prev->front, next->rear);
	}
	static void setNeighbors(MovingParticle* p, MovingParticle* prev, MovingParticle* next, MovingFront rear, MovingFront front)
	{
		p->prev = prev;
		p->next = next;
		p->prev->next = p;
		p->next->prev = p;
		p->rear = rear;
		p->front = front;
		/*p->reflexive = GetVisualAngle2(p->p.m_X - rear.dir.x, p->p.m_Y - p->rear.dir.y,
			p->p.m_X + p->front.dir.x, p->p.m_Y + p->front.dir.y,
			p->p.m_X, p->p.m_Y);*/
	}
	static vector<MovingParticle*> vectorize(MovingParticle* p);
	static vector<MovingParticle*> extractPath(MovingParticle* p, MovingParticle* q);
	static bool updatePolygon(MovingParticle* p);
	static MovingParticle* getNextEvent();
	static vector<MovingParticle*> getNextEvents(float eps = 1.0e-6);
	//static vector<vector<CParticleF>> takeSnapshots();
	//static void removeUnstable(); //remove kinks and duplicates.
	static void quickFinish(); //remove polygons with less than 4 vertices.
	static bool sanityCheck();
	static vector<vector<MovingParticle*>> clusterParticles();
	//static void _traceBack(MovingParticle* p, vector<MovingParticle*>& trace, set<MovingParticle*>& pset);
	static vector<MovingParticle*> MovingParticle::traceBackPolygon(vector<MovingParticle*>& particles);
	static vector<vector<MovingParticle*>> closedRegions(vector<MovingParticle*>& points);
	static vector<vector<MovingParticle*>> closedRegions2(vector<MovingParticle*>& points);
	static pair<MovingParticle*, float> findIntersection(MovingParticle* p, MovingParticle* q);
	static float intersectSideAndLine(const MovingParticle* p, const MovingParticle* q, const MovingParticle* r);
	static MovingParticle* traceAndHandleUnstable(MovingParticle* p, MovingParticle* parent);
	static bool correctOvershoot(MovingParticle* p, MovingParticle* q, pair<float, float> param);

private:
	MovingParticle(StationaryParticle* p = NULL, MovingParticleType t = Unknown, float tm = 0.0f)
	{
		init_particle = p;
		init_particle->derived.insert(this);
		this->p = p->getP();
		p0 = p->getP();
		v[0] = v[1] = std::numeric_limits<float>::quiet_NaN();
		next = NULL;
		prev = NULL;
		parents[0] = parents[1] = parents[2] = parents[3] = NULL;
		//children[0] = children[1] = NULL;
		id = _id++;
		type = t;
		created = tm;
		time = created;
		bActive = true;
		bInitialized = false;
		bUnstable = false;
		//bNeedUpdate = true;
		bLeaf = false;
		init_event_time = -1;
		reflexive = std::numeric_limits<float>::quiet_NaN();
		event = EventStruct(std::numeric_limits<float>::infinity(), UnknownEvent, this);
		polygon = NULL;
	}

	float MovingParticle::_splitTime(const MovingParticle* q, const MovingParticle* r, float eps = 1.0e-3) const;
	bool _onSideAt(const MovingParticle* q, float t, float eps = 1.0e-3) const;
	void _setParents(EventStruct cause, bool bSide);

	StationaryParticle* init_particle;
	CParticleF p0;
	CParticleF p;
	MovingParticle* next;
	MovingParticle* prev;
	MovingParticle* parents[NumParents];
	//MovingParticle* children[2];
	MovingFront rear;
	MovingFront front;

	int id;
	MovingParticleType type;
	float created; //time this is created.
	float time;
	float init_event_time;
	bool bInitialized; //false until its velocity is set.
	bool bActive; //true if it is still moving.
	bool bUnstable; //true if the bisector angle was too tight for accurate velocity calculation.
	bool bLeaf; //true if this is a leaf of the tree
	//bool bNeedUpdate; //true if the event may need to be udpated.
	float reflexive; //measure of reflexivenessbReflexive
	EventStruct event;
	float v[2];
	Polygon* polygon;
	static int _id;
};

struct ParticleFactory
{
	static ParticleFactory& getInstance()
	{
		static ParticleFactory* _instance = NULL;
		if (_instance == NULL)
		{
			_instance = new ParticleFactory();
		}
		return *_instance;

	}
	MovingParticle* makeParticle(StationaryParticle* p, MovingParticleType type, float tm)
	{
		MovingParticle* particle = new MovingParticle(p, type, tm);
		particles.push_back(particle);
		activeSet.insert(particle);
		pmap[particle->id] = particle;
		if (particle->id == 51)
		{
			particle->id += 0;
		}
		return particle;
	}
	bool inactivate(MovingParticle* p)
	{
		set<MovingParticle*>::iterator it = activeSet.find(p);
		if (it == activeSet.end())
		{
			return false;
		}
		else
		{

			(*it)->bActive = false;
			activeSet.erase(it);
			return true;
		}
	}
	//when initialized from snapshots, some future particles have to be physically removed from memory.
	bool remove(MovingParticle* p)
	{
		particles.erase(find(particles.begin(), particles.end(), p));
		pmap.erase(p->id);
		set<MovingParticle*>::iterator it = find(activeSet.begin(), activeSet.end(), p);
		if (it != activeSet.end())
		{
			activeSet.erase(it);
		}
		delete p;
		return true;
	}
	void clean()
	{
		activeSet.clear();
		for (int i = 0; i<particles.size(); ++i)
		{
			delete particles[i];
		}
		particles.clear();
		pmap.clear();
		MovingParticle::_id = 0;
	}
	MovingParticle* get(int id)
	{
		return pmap[id];
	}
	MovingParticle* getNext()
	{
		float t = numeric_limits<float>::infinity();
		MovingParticle* p = NULL;
		for (set<MovingParticle*>::iterator it = activeSet.begin(); it != activeSet.end(); ++it)
		{
			MovingParticle* q = *it;
			if (q->event.t < t)
			{
				t = q->event.t;
				p = q;
			}
		}
		return p;
	}
	vector<MovingParticle*> particles;
	set<MovingParticle*> activeSet;
	map<int, MovingParticle*> pmap;
	set<MovingParticle*> updateQueue;
private:
	ParticleFactory()
	{
		MovingParticle::_id = 0;
	}
	ParticleFactory(ParticleFactory& f)
	{
	}
	ParticleFactory operator=(ParticleFactory& f)
	{
	}
	~ParticleFactory()
	{
		clean();
	}
};

/*
Utility functions
*/

bool
calculateBisectorVelocity(CParticleF a, CParticleF o, CParticleF b, float& vx, float& vy);

vector<MovingParticle*> vectorize(MovingParticle* p);

/*
true if the vector of points are in clock-wise direction
*/
int clockWise(vector<MovingParticle*>& particles);

