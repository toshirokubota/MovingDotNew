#include <MovingParticle.h>
#include <vector>
#include <set>
#include <limits>
using namespace std;
#include <mex.h> //need for premature exit when something goes wrong.

#include <szmexutilitytemplate.h>
#include <szMiscOperations.h>
#include <IntersectionConvexPolygons.h>
#include <Graph.h>
#include <ParticleDirection.h>
#include <Polygon.h>

const float MP_EPSILON = 1.0e-2;
const float ENDPOINT_SPEED = 1.0f; //Make it >1.0 to let the end points move faster.

//enum MovingParticleType { Unknown, Initial, Regular, Merge, Collide, Split1, Split2, Collide1, Collide2, Axis, Dummy };
MovingParticleType int2ParticleType(int i)
{
	MovingParticleType type = Unknown;
	switch (i)
	{
	case 1:
		type = Initial;
		break;
	case 2:
		type = Regular;
		break;
	case 3:
		type = Merge;
		break;
	case 4:
		type = Split;
		break;
	case 5:
		type = Collide;
		break;
	case 6:
		type = Axis;
		break;
	case 7:
		type = Dummy;
		break;
	}
	return type;
}

void
MovingParticle::print(char* tab) const
{
	printf("%s%d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %d\n", tab, id, p0.m_X, p0.m_Y, p.m_X, p.m_Y, v[0], v[1], init_particle->getId());
}

void
MovingParticle::printParentTree(char* tab) const
{
	//if (this->parents[0] == NULL && this->parents[1] == NULL)
	{
		printf("%s%d %3.3f %3.3f %3.3f %3.3f %d\n", tab, id, p0.m_X, p0.m_Y, p.m_X, p.m_Y, type);
	}
	//else
	{
		for (int i = 0; i < MovingParticle::NumParents; ++i)
		{
			if (this->parents[i] != NULL)
			{
				char tab2[256];
				strcpy(tab2, tab);
				strcat(tab2, "\t");
				this->parents[i]->printParentTree(tab2);
			}
		}
	}
}

vector<float>
MovingParticle::dump2vector()
{
	vector<float> v;
	v.push_back(id);
	v.push_back(p0.m_X);
	v.push_back(p0.m_Y);
	v.push_back(p.m_X);
	v.push_back(p.m_Y);
	v.push_back(this->v[0]);
	v.push_back(this->v[1]);
	v.push_back(prev == NULL ? -1 : prev->id);
	v.push_back(next == NULL ? -1 : next->id);
	v.push_back(type);
	v.push_back(created);
	v.push_back(time);
	v.push_back(reflexive);
	v.push_back(bActive ? 1.0f : 0.0f);
	v.push_back(bUnstable ? 1.0f : 0.0f);
	v.push_back(parents[0] == NULL ? -1 : parents[0]->id);
	v.push_back(parents[1] == NULL ? -1 : parents[1]->id);
	v.push_back(event.type);
	v.push_back(event.q == NULL ? -1.0f : event.q->id);
	v.push_back(event.r == NULL ? -1.0f : event.r->id);
	v.push_back(event.t);
	v.push_back(front.dir.x);
	v.push_back(front.dir.y);
	v.push_back(front.speed);
	v.push_back(rear.dir.x);
	v.push_back(rear.dir.y);
	v.push_back(rear.speed);
	return v;
}

/*
A reflex particle is non convex one that can split a side of a polygon.
*/
bool 
MovingParticle::isReflex() //check if it is concave (reflexive) that allows splitting of a side;
{
	//return GetVisualAngle2(prev->p.m_X, prev->p.m_Y, next->p.m_X, next->p.m_Y, p.m_X, p.m_Y) <= 0;
	//return this->reflexive <= 0.0f; ///TK: should 0 be included? 8/28/2015
	if (reflexive != reflexive)
	{
		reflexive = GetVisualAngle2(p.m_X - rear.dir.x, p.m_Y - rear.dir.y,
			p.m_X + front.dir.x, p.m_Y + front.dir.y,
			p.m_X, p.m_Y);
	}
	return reflexive <= 0.0f;
}

//find the angle of propagating front p-q.
float 
MovingParticle::frontPropAngle(MovingParticle* p, MovingParticle* q)
{
	CParticleF p0 = p->move(1);
	CParticleF p1 = p->p;
	CParticleF p2 = q->p;
	CParticleF p3 = q->move(1);

	return GetVisualAngle2(p0.m_X, p0.m_Y, p2.m_X, p2.m_Y, p1.m_X, p1.m_Y) + 
		GetVisualAngle2(p1.m_X, p1.m_Y, p3.m_X, p3.m_Y, p2.m_X, p2.m_Y);

}

vector<MovingParticle*>
MovingParticle::vectorize(MovingParticle* start)
{
	/*vector<MovingParticle*> tr;
	MovingParticle* p = start;
	set<MovingParticle*> pset;
	bool success = true;
	do
	{
		tr.push_back(p);
		if (pset.find(p) != pset.end()) //premature loop is found
		{
			success = false;
			break;
		}
		if (p->next == NULL)
		{
			success = false;
			break;
		}
		pset.insert(p);
		p = p->next;
	} while (p != start);
	if (success == false)
	{
		printf("failed vectorization data.\n");
		for (int i = 0; i<tr.size(); ++i)
		{
			printf("%d %f %f %d\n", i + 1, tr[i]->p.m_X, tr[i]->p.m_Y, tr[i]->id);
		}
		mexErrMsgTxt("vectorize: failed to vectorize a shape.");
		//tr.clear();
	}
	return tr;*/
	return extractPath(start, start);
}

vector<MovingParticle*>
MovingParticle::extractPath(MovingParticle* start, MovingParticle* end)
{
	vector<MovingParticle*> tr;
	MovingParticle* p = start;
	set<MovingParticle*> pset;
	bool success = true;
	do
	{
		tr.push_back(p);
		if (pset.find(p) != pset.end()) //premature loop is found
		{
			success = false;
			break;
		}
		if (p->next == NULL)
		{
			success = false;
			break;
		}
		pset.insert(p);
		p = p->next;
	} while (p != end);
	if (success == false)
	{
		printf("failed vectorization data.\n");
		for (int i = 0; i<tr.size(); ++i)
		{
			printf("%d %f %f %d\n", i + 1, tr[i]->p.m_X, tr[i]->p.m_Y, tr[i]->id);
		}
		mexErrMsgTxt("vectorize: failed to vectorize a shape.");
		//tr.clear();
	}
	return tr;
}

bool
MovingParticle::updatePolygon(MovingParticle* p)
{
	vector<MovingParticle*> vp = vectorize(p);
	PolygonFactory& factory = PolygonFactory::getInstance();
	Polygon* polygon = factory.makePolygon(vp, p->time);
	for (int i = 0; i < vp.size(); ++i)
	{
		vp[i]->setPolygon(polygon);
	}
	return true;
}

//find the particle with the next event in line.
MovingParticle*
MovingParticle::getNextEvent()
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	MovingParticle* p = NULL;
	float time = std::numeric_limits<float>::infinity();
	for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); it++)
	{
		EventStruct ev = (*it)->event;
		if (ev.t < time || (p != NULL && ev.t == time && ev.p->id > p->id))
		{
			p = *it;
			time = ev.t;
		}
	}
	return p;
}

/*
Takes care of end-points (leaf nodes) and set rear and front directions.
This should be called at the very beginning of the simulation after forrests are traced.
*/
bool
MovingParticle::initializeVelocity()
{
	if (bInitialized) return false; //cannot be set more than once.

	float eps = 1.0e-6;
	int mode = 0; //0: normal, 1: prev-leaf, 2: next-leaf
	if (Distance(p, prev->p) < eps && Distance(p, next->p) < eps) //a single point. Cann't continue
	{
		v[0] = v[1] = 0;
		return false;
	}
	else
	{
		if (Distance(p, next->p) > eps)
		{
			CParticleF d = NormalizedDirection(next->p, p);
			front.dir.x = d.m_X;
			front.dir.y = d.m_Y;
			front.speed = 1.0f;
		}
		else //p==next
		{
			//leaf
			mode = 2;
			float ang = GetVisualDirection(p.m_X, p.m_Y, prev->p.m_X, prev->p.m_Y) - PI / 2.0;
			CParticleF d = perpDirection(prev->p, p);
			front.dir.x = cos(ang);
			front.dir.y = sin(ang);
			front.speed = ENDPOINT_SPEED;
			bLeaf = true;
		}
		if (Distance(p, prev->p) > eps)
		{
			CParticleF d = NormalizedDirection(p, prev->p);
			rear.dir.x = d.m_X;
			rear.dir.y = d.m_Y;
			rear.speed = 1.0f;
		}
		else//p==prev
		{
			//leaf
			mode = 1;
			float ang = GetVisualDirection(p.m_X, p.m_Y, next->p.m_X, next->p.m_Y) + PI / 2.0;
			CParticleF d = perpDirection(p, next->p);
			rear.dir.x = -cos(ang);
			rear.dir.y = -sin(ang);
			rear.speed = ENDPOINT_SPEED;
			bLeaf = true;
		}
		bool b = calculateVelocityR();
		return b;
	}
}

bool
MovingParticle::calculateVelocityR()
{
	if (bInitialized) return false; //cannot be set more than once.
	if (next == prev) return false; //kink
	bInitialized = true;

	CParticleF a1(p.m_X - rear.speed*rear.dir.y, p.m_Y + rear.speed*rear.dir.x);
	CParticleF a2(a1.m_X + rear.dir.x, a1.m_Y + rear.dir.y);
	CParticleF b1(p.m_X - front.speed*front.dir.y, p.m_Y + front.speed*front.dir.x);
	CParticleF b2(b1.m_X + front.dir.x, b1.m_Y + front.dir.y);
	pair<float, float> pr = _IntersectConvexPolygon::intersect(a1, a2, b1, b2);
	if (pr.first != pr.first)
	{
		//on top of a straight line segment
		ParticleDirection pd = rear.normal();
		float speed = (rear.speed + front.speed) / 2.0;
		v[0] = speed * pd.x;
		v[1] = speed * pd.y;
		return true;
	}
	else if (Abs(pr.first) > 1000.0)
	{
		bUnstable = true;
		return false;
	}
	else
	{
		CParticleF c((1.0 - pr.first) * a1.m_X + pr.first * a2.m_X, (1.0 - pr.first) * a1.m_Y + pr.first * a2.m_Y);
		v[0] = c.m_X - p.m_X;
		v[1] = c.m_Y - p.m_Y;
		return true;
	}
}

/*bool
MovingParticle::calculateVelocity(bool bleaf)
{
	CParticleF b(p.m_X + front.x, p.m_Y + front.y);
	CParticleF a(p.m_X - rear.x, p.m_Y - rear.y);
	CParticleF bs = bisector(p, a, b);
	if (bleaf)
	{
		if (Distance(p0, next->p0) <= 1.0e-3)
		{
			b = CParticleF(p.m_X - bs.m_X, p.m_Y - bs.m_Y);
			bs = bisector(p, a, b);
		}
		else if (Distance(p0, prev->p0) <= 1.0e-3)
		{
			a = CParticleF(p.m_X - bs.m_X, p.m_Y - bs.m_Y);
			bs = bisector(p, a, b);
		}
	}
	double ang = GetVisualAngle2(a.m_X, a.m_Y, b.m_X, b.m_Y, p.m_X, p.m_Y);
	double cs = cos((PI - Abs(ang)) / 2.0);
	bool bret = true;
	float eps = 0.0001;
	if (cs < eps)
	{
		cs = eps;
		bret = false;
	}
	double len = 1.0f / cs;
	bs.m_X *= len;
	bs.m_Y *= len;
	v[0] = (float)bs.m_X; 
	v[1] = (float)bs.m_Y;
	return bret;
}*/

int clockWise(vector<MovingParticle*>& particles)
{
	vector<CParticleF> pnts;
	for (int i = 0; i<particles.size(); ++i)
	{
		pnts.push_back(particles[i]->getP());
	}
	return ClockWise(pnts);
}

/*
This utility method check when a moving particle (P) with a moving line of Q and its next.
It returns the time of the split/collision.
*/
float 
MovingParticle::_splitTime(const MovingParticle* q, const MovingParticle* r, float eps) const
{
	ParticleDirection u = q->front.normal();
	float dq = q->v[0] * u.x + q->v[1] * u.y;
	if (dq < 0)
	{
		u.x = -u.x;
		u.y = -u.y;
		dq = -dq;
	}
	float dp = (u.x * v[0] + u.y * v[1]);
	if (dp >= 0) return std::numeric_limits<float>::infinity(); //moving along the same direction

	ParticleDirection ve((1 - dp)*u.x, (1 - dp)*u.y);

	CParticleF y = Closest2Line(q->p, CParticleF(q->p.m_X + q->front.dir.x, q->p.m_Y + q->front.dir.y), p);
	float deriv = (y.m_X - p.m_X)*ve.x + (y.m_Y - p.m_Y) * ve.y;
	if (deriv >= 0) return std::numeric_limits<float>::infinity(); //moving away

	float dval = Distance(p, y);
	return dval / (dq - dp);
}

/*
*/
/*
Check if the particle P will be located on the side formed by Q and its next at time t from the current time.
*/
bool
MovingParticle::_onSideAt(const MovingParticle* q, float t, float eps) const
{
	CParticleF p2 = this->move(t);
	CParticleF q2 = q->move(t);
	CParticleF r2 = q->next->move(t);
	float d = Distance2LineSegment(q2, r2, p2);
	return d <= MP_EPSILON; //TK NEED TO FIX THIS.
	//pair<float, float> param = _IntersectConvexPolygon::intersect(p, p2, q2, r2);
	//return param.second >= 0.0f && param.second <= 1.0f;
}

/*
For a particle created by an event (CAUSE), set its parents. 
If bSide is false, the natural order is reversed.
*/
void
MovingParticle::_setParents(EventStruct cause, bool bSide)
{
	MovingParticle* pe = (MovingParticle*)cause.p;
	MovingParticle* qe = (MovingParticle*)cause.q;
	MovingParticle* re = (MovingParticle*)cause.r;
	if (cause.type == EdgeEvent)
	{
		this->parents[0] = pe;
		this->parents[1] = qe;
	}
	else if (cause.type == SplitEvent)
	{
		float dpq = Distance(p, qe->p);
		float dpr = Distance(p, re->p);
		MovingParticle* x = NULL;
		if (dpq <= dpr)
		{
			x = qe;
		}
		else
		{
			x = re;
		}
		if (bSide)
		{
			this->parents[0] = pe;
			this->parents[1] = x;
		}
		else
		{
			this->parents[0] = x;
			this->parents[1] = pe;
		}
	}
	else if (cause.type == CollisionEvent)
	{
		if (bSide)
		{
			this->parents[0] = pe;
			this->parents[1] = qe;
		}
		else 
		{
			this->parents[0] = qe;
			this->parents[1] = pe;
		}
	}
}

/*
Find the time to hit for the 1st particle (P) to a slanted plane formed with Q and R.
Returns Inf if no intersection within the bounded polygon.
*/
float
MovingParticle::intersectSideAndLine(const MovingParticle* p, const MovingParticle* q, const MovingParticle* r)
{
	float t = std::numeric_limits<float>::infinity();
	float t0 = p->_splitTime(q, r);
	if (t0 > 0)
	{
		if (p->_onSideAt(q, t0))
		{
			t = t0;
		}
	}
	return t;
}

EventStruct
MovingParticle::findNextEdgeEvent() const
{
	CParticleF p2(p.m_X + v[0], p.m_Y + v[1]);
	CParticleF q0 = next->p;
	CParticleF q2(q0.m_X + next->v[0], q0.m_Y + next->v[1]);
	pair<float,float> param = _IntersectConvexPolygon::intersect(p, p2, q0, q2);
	EventStruct ev(std::numeric_limits<float>::infinity(), EdgeEvent, this, this->next);
	if (param.first >= 0.0f && param.second >= 0.0f)
	{
		ev.t = time + param.first;
	}
	/*if (id == 338)
	{
		printf("%d: %f %f %f %f %f %f %f %f %f %f\n", id, param.first, param.second, p.m_X, p.m_Y, p2.m_X, p2.m_Y, q0.m_X, q0.m_Y, q2.m_X, q2.m_Y);
	}*/
	return ev;
}

EventStruct
MovingParticle::findNextSplitEvent() const
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	EventStruct ev(std::numeric_limits<float>::infinity(), SplitEvent, this);
	for (set<MovingParticle*>::iterator j = factory.activeSet.begin(); j != factory.activeSet.end(); ++j)
	{
		const MovingParticle* q = *j; 
		const MovingParticle* r = q->next;
		if (id == 373 && (q->id == 99))
			ev.t += 0;
		if (this == q || this->next == q || this==r || this->prev==r) continue;
		/*if ((Abs(this->created - q->created) < 1.0e-8 && Distance(this->p, q->p) < 1.0e-3) ||
			(Abs(this->created - r->created) < 1.0e-8 && Distance(this->p, r->p) < 1.0e-3))
		{
			//Nov. 15th, 2014
			//two particles created at the same time at the same place.
			//they cannot be splitting each other. However, because of numerical precision, we may get a very small positive 
			//collision time.
			continue;
		}*/

		if (init_particle == q->init_particle || init_particle == r->init_particle) continue; //they are technically the same particle.

		//if (eta <= 0) continue;
		float t0 = intersectSideAndLine(this, q, r); // period - p->time);
		if (t0 >= 0)
		{
			float t = t0 + time;
			if (t < ev.t)
			{
				ev.t = t;
				ev.q = q;
				ev.r = r;
			}
		}
	}
	//if the event location is close to a vertex of being split, then make it as Collision
	if(ev.q != NULL)
	{
		float t0 = ev.t - time;
		CParticleF p0 = move(t0);
		float dq = Distance(p0, ev.q->move(t0));
		float dr = Distance(p0, ev.r->move(t0));
		if ( dq < MP_EPSILON || dr < MP_EPSILON)
		{
			ev.type = CollisionEvent;
			if (dr < dq)
			{
				ev.q = ev.r;
				ev.r = ev.q->next;
			}
		}
	}

	return ev;
}

/*
restrict serch of the next event within the given set. The set typically is a set of particles in the current polygon.
*/
EventStruct
MovingParticle::findNextSplitEvent(const vector<MovingParticle*>& vp) const
{
	EventStruct ev(std::numeric_limits<float>::infinity(), SplitEvent, this);
	for (int j = 0; j<vp.size(); ++j)
	{
		const MovingParticle* q = vp[j];
		if (id == 532 && q->id == 487)
			j += 0;
		const MovingParticle* r = q->next;
		if (this == q || this->next == q || this == r || this->prev == r) continue;
		if (init_particle == q->init_particle || init_particle == r->init_particle) continue; //they are technically the same particle.

		float t0 = intersectSideAndLine(this, q, r); // period - p->time);
		if (t0 >= 0)
		{
			float t = t0 + time;
			if (t < ev.t)
			{
				ev.t = t;
				ev.q = q;
				ev.r = r;
			}
		}
	}
	//if the event location is close to a vertex of being split, then make it as Collision
	if (ev.q != NULL)
	{
		float t0 = ev.t - time;
		CParticleF p0 = move(t0);
		float dq = Distance(p0, ev.q->move(t0));
		float dr = Distance(p0, ev.r->move(t0));
		if (dq < MP_EPSILON || dr < MP_EPSILON)
		{
			ev.type = CollisionEvent;
			if (dr < dq)
			{
				ev.q = ev.r;
				ev.r = ev.q->next;
			}
		}
	}

	return ev;
}

bool
MovingParticle::updateEvent()
{
	//check if the event needs to be updated. If not, then simply return.
	if (this->bUnstable) return false; //unstable particle.
	if (this->bActive == false) return  false;
	bool bFirst = event.type == UnknownEvent;
	if (id == 255 && next->id==765 || id==765 && next->id==187)
		id += 0;

	if (event.t >= std::numeric_limits<float>::infinity())
	{
		EventStruct ev2;
		if (isReflex()) 
		{
			ev2 = findNextSplitEvent();
		}
		EventStruct ev1 = findNextEdgeEvent();
		if (ev1.t <= ev2.t) 
		{
			event = ev1;
		}
		else
		{
			event = ev2;
		}
	}
	else 
	{
		EventStruct ev1;
		EventStruct ev2;
		EventStruct ev3;
		float t0 = event.t;
		bool bChanged = false;
		if (isReflex() && (event.q->next != event.r || event.q->isActive() == false || event.r->isActive() == false))
		{
			ev2 = findNextSplitEvent();
			ev1 = findNextEdgeEvent();
			bChanged = true;
		}
		if (event.type == EdgeEvent && (this->next != event.q || event.q->isActive() == false))
		{
			ev1 = findNextEdgeEvent();
			bChanged = true;
		}
		if (this->next->p == this->next->p0) //newly created neighbor
		{
			ev3 = findNextEdgeEvent();
			if (ev3.t < t0)
			{
				bChanged = true;
			}
		}
		if (bChanged)
		{
			if (ev1.t < ev2.t && ev1.t < ev3.t) // && ev1.t < std::numeric_limits<float>::infinity())
			{
				event = ev1;
			}
			else if (ev2.t < ev1.t && ev2.t < ev3.t)
			{
				event = ev2;
			}
			else
			{
				event = ev3;
			}
		}
		//event.t = Min(t0, event.t); //event time cannot increase.
	}

	if (bFirst) {
		this->init_event_time = event.t;
	}

	return true;
}


pair<MovingParticle*, MovingParticle*>
MovingParticle::applyEvent()
{
	float eps = 1.0e-3;
	MovingParticle* p = this;
	MovingParticle* q = (MovingParticle*)event.q;
	MovingParticle* r = (MovingParticle*)event.r;
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	StationaryParticle* sp = sfactory.makeParticle(p->p);
	MovingParticle* pnew[2] = { NULL, NULL };
	if (event.type == EdgeEvent)
	{
		pnew[0] = factory.makeParticle(sp, Collide, event.t);
		setNeighbors(pnew[0], p->prev, p->next->next);
		factory.inactivate(p);
		factory.inactivate(q);
	}
	else if (event.type == SplitEvent)
	{
		pnew[0] = factory.makeParticle(sp, Split, event.t);
		pnew[1] = factory.makeParticle(sp, Split, event.t);

		setNeighbors(pnew[0], event.p->prev, r);
		setNeighbors(pnew[1], q, event.p->next);
		factory.inactivate(p);
	}
	else if (event.type == CollisionEvent) //colliding at two vertices
	{
		pnew[0] = factory.makeParticle(sp, Collide, event.t);
		pnew[1] = factory.makeParticle(sp, Collide, event.t);
		MovingParticle* x = (MovingParticle*)event.q;
		MovingParticle* xn = x->next;
		MovingParticle* xp = x->prev;
		setNeighbors(pnew[0], event.p->prev, xn);
		setNeighbors(pnew[1], xp, event.p->next);
		factory.inactivate(p);
		factory.inactivate(x);
	}

	for (int i = 0; i < 2; ++i)
	{
		if (pnew[i] == NULL) continue;

		pnew[i]->_setParents(event, i==0); //set parents of the new particle
		if (pnew[i]->calculateVelocityR() == false)
		{
			pnew[i]->bUnstable = true;
			/*pnew[i]->parents[3] = pnew[i]->parents[1];
			pnew[i]->parents[1] = pnew[i]->prev;
			pnew[i]->parents[2] = pnew[i]->next;*/
			MovingParticle* ptmp = pnew[i];
			pnew[i] = traceAndHandleUnstable(pnew[i], pnew[i==0? 1: 0]);
			factory.inactivate(ptmp);
		}
		updatePolygon(pnew[i]);
	}

	return 	pair<MovingParticle*, MovingParticle*> (pnew[0], pnew[1]);
}

MovingParticle*
MovingParticle::traceAndHandleUnstable(MovingParticle* p, MovingParticle* parent)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	if (p->bActive == false)
	{
		return p;
	}
	if (p->prev == p->next)
	{
		/*MovingParticle* pnew = factory.makeParticle(sfactory.makeParticle(p->p), Dummy, p->time);
		pnew->setNeighbors(pnew, pnew, pnew); //pure dummy
		pnew->parents[0] = p;
		pnew->parents[1] = parent;
		factory.inactivate(pnew);
		return pnew;*/
		//return p; //there are only two particles in this chain. quickFinish will take care of this.
	}
	float dp = Distance(p->p, p->prev->p);
	float dn = Distance(p->p, p->next->p);
	MovingParticle* q = NULL;
	if (dp < dn)
	{
		q = factory.makeParticle(sfactory.makeParticle(p->prev->p), Split, p->time);
		q->setNeighbors(q, p->prev->prev, p->next);
		q->parents[0] = p->parents[0];
		q->parents[1] = p->parents[1];
		factory.inactivate(p->prev);
	}
	else
	{
		q = factory.makeParticle(sfactory.makeParticle(p->next->p), Split, p->time);
		q->setNeighbors(q, p->prev, p->next->next);
		q->parents[0] = p->parents[0];
		q->parents[1] = p->parents[1];
		factory.inactivate(p->next);
	}

	if (q->calculateVelocityR() == false)
	{
		//MovingParticle* pnew = factory.makeParticle(sfactory.makeParticle(q->p), Axis, q->time);
		MovingParticle* q2 = traceAndHandleUnstable(q, p);
		factory.inactivate(q);
		//factory.inactivate(pnew);
		return q2;
	}
	else {
		return q;
	}
}

void
MovingParticle::quickFinish()
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	set<MovingParticle*> pset;
	while (true)
	{
		bool bdone = true;
		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (pset.find(p) == pset.end())
			{
				vector<MovingParticle*> vp = vectorize(p);
				for (int i = 0; i < vp.size(); ++i)
				{
					pset.insert(vp[i]);
				}
				if (vp.size() == 3)
				{
					CParticleF x;
					for (int i = 0; i < vp.size(); ++i)
					{
						x.m_X += vp[i]->p.m_X/3;
						x.m_Y += vp[i]->p.m_Y/3;
					}
					MovingParticle* pnew[2];
					pnew[0] = factory.makeParticle(sfactory.makeParticle(x), Dummy, vp[0]->time);
					pnew[1] = factory.makeParticle(sfactory.makeParticle(x), Axis, vp[0]->time);
					pnew[0]->parents[0] = vp[0];
					pnew[0]->parents[1] = vp[1];
					pnew[1]->parents[0] = pnew[0];
					pnew[1]->parents[1] = vp[2];
					for (int i = 0; i < 2; ++i)
					{
						factory.inactivate(pnew[i]);
					}
				}
				if (vp.size() <= 3)
				{
					for (int i = 0; i < vp.size(); ++i)
					{
						factory.inactivate(vp[i]);
					}
				}
				bdone = false;
				break; //active set has changed. the interator needs to be initialized again.
			}
		}
		if (bdone) break;
	}
}

/*
find a segment in a polygon wich Q where [p, p->next] intersects.
It returns the segment [q, q->next] and a parameter t, such that the intersection point is
(1-t)*q + t*q->next.
*/
pair<MovingParticle*, float>
MovingParticle::findIntersection(MovingParticle* p, MovingParticle* q)
{
	pair<MovingParticle*, float> info(NULL, 0.0f);
	vector<MovingParticle*> vp = vectorize(q);
	for (int i = 0; i<vp.size(); ++i)
	{
		MovingParticle* r = vp[i];
		pair<float, float> param = _IntersectConvexPolygon::intersect(p->p, p->next->p, r->p, r->next->p);
		if (param.first > 0 && param.first < 1.0f &&param.second >= 0 && param.second <= 1.0f)
		{
			info.first = r;
			info.second = param.second;
			break;
		}
	} 
	return info;
}

bool 
MovingParticle::correctOvershoot(MovingParticle* p, MovingParticle* q, pair<float, float> param)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	float time = p->time;
	if (param.first > 0.0f && param.first<1.0f && param.second>0.0f && param.second < 1.0f)
	{
		if (p->next->next == q) //collision
		{
			MovingParticle* p2 = p->next;
			float t = param.first;
			CParticleF p0(p->p.m_X*(1.0 - t) + p2->p.m_X*t, p->p.m_Y*(1.0 - t) + p2->p.m_Y*t);
			MovingParticle* y = factory.makeParticle(sfactory.makeParticle(p0), Collide, time);
			setNeighbors(y, p, q->next);
			y->calculateVelocityR();
			y->parents[0] = p2; //TK?
			y->parents[1] = q;  //TK?
			factory.inactivate(p2);
			factory.inactivate(q);
			vector<MovingParticle*> vp1 = vectorize(y);
		}
		else //split
		{
			MovingParticle* rs[2];
			if (param.first < 0.5) //intersection closer to p than p->next
			{
				rs[0] = p->prev;
				rs[1] = p;
			}
			else
			{
				rs[0] = p;
				rs[1] = p->next;
			}
			pair<MovingParticle*, float> isec[2]; //intersection information
			for (int k = 0; k < 2; ++k)
			{
				isec[k] = findIntersection(rs[k], q);
			}
			if (isec[0].first == NULL || isec[1].first == NULL)
			{
				printf("correctOvershoot: %d <=> %d, %d <=> %d\n",
					rs[0]->id, (isec[0].first == NULL ? -1 : isec[0].first->id),
					rs[1]->id, (isec[1].first == NULL ? -1 : isec[1].first->id));
				return false;
			}
			MovingParticle* pnew[2];
			for (int k = 0; k < 2; ++k)
			{
				float t = isec[k].second;
				MovingParticle* x = isec[k].first;
				CParticleF p0(x->p.m_X*(1.0 - t) + x->next->p.m_X*t, x->p.m_Y*(1.0 - t) + x->next->p.m_Y*t);
				pnew[k] = factory.makeParticle(sfactory.makeParticle(p0), Split, time);
			}

			//inactivate some before changing the neighborhood configuration
			factory.inactivate(rs[1]);
			MovingParticle* tmp = isec[1].first;
			while (tmp != isec[0].first)
			{
				factory.inactivate(tmp->next);
				tmp = tmp->next;
			}
			//set neighborhoods, velocity and parents.
			setNeighbors(pnew[0], rs[0], isec[0].first->next);
			setNeighbors(pnew[1], isec[1].first, rs[1]->next);
			vector<MovingParticle*> vp1 = vectorize(pnew[0]);
			vector<MovingParticle*> vp2 = vectorize(pnew[1]);
			pnew[0]->parents[0] = rs[1]; //TK?
			pnew[1]->parents[1] = rs[1]; //TK?
			for (int k = 0; k < 2; ++k)
			{
				pnew[k]->calculateVelocityR();
			}
		}
	}
	return true;
}

bool
MovingParticle::sanityCheck()
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	while (true)
	{
		bool bdone = true;
		for (int i = 0; i < factory.particles.size() && bdone; ++i)
		{
			MovingParticle* p = factory.particles[i];
			if (p->isActive() == false) continue;
			MovingParticle* p2 = p->next;

			for (int j = 0; j < factory.particles.size(); ++j)
			{
				MovingParticle* q = factory.particles[j];
				if (q->isActive() == false) continue;
				MovingParticle* q2 = q->next;
				if (p2 == q || q2 == p) continue;
				pair<float, float> param = _IntersectConvexPolygon::intersect(p->p, p2->p, q->p, q2->p);
				if (param.first > 0.0f && param.first<1.0f && param.second>0.0f && param.second < 1.0f)
				{
					printf("Correction: %d-%d vs %d-%d (%f, %f)\n", p->id, p2->id, q->id, q2->id, param.first, param.second);
					bool b = correctOvershoot(p, q, param);
					if (b == false)
					{
						printf("Correction failed.\n");
						//return false;
					}
					else
					{
						bdone = false;
						break;
					}
				}
			}
		}

		if (bdone) break;
	}
	return true;
}

vector<vector<MovingParticle*>> 
MovingParticle::clusterParticles()
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<vector<MovingParticle*>> clusters;
	set<MovingParticle*> pset;
	while (true)
	{
		bool bdone = true;
		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (pset.find(p) == pset.end())
			{
				vector<MovingParticle*> vp = MovingParticle::vectorize(p);
				for (int i = 0; i < vp.size(); ++i)
				{
					pset.insert(vp[i]);
				}
				bdone = false;
				clusters.push_back(vp);
			}
		}
		if (bdone) break;
	}
	return clusters;
}

void
_traceBack(MovingParticle* p, vector<MovingParticle*>& trace, set<MovingParticle*>& pset)
{
	if (p == NULL){
		p = NULL;
	}
	pset.insert(p);
	if (p->getType() == Initial)
	{
		trace.push_back(p);
	}
	else
	{
		for (int i = 0; i < MovingParticle::NumParents; ++i)
		{
			if (p->getParent(i) != NULL)
			{
				_traceBack(p->getParent(i), trace, pset);
			}
		}
	}
}

/*
Given a vectorized particles, it returns a series of initial particles by tracing particles backward in time.
*/
vector<MovingParticle*>
MovingParticle::traceBackPolygon(vector<MovingParticle*>& particles)
{
	vector<MovingParticle*> trace;
	set<MovingParticle*> pset;
	for (int i = 0; i<particles.size(); ++i)
	{
		if (pset.find(particles[i]) == pset.end())
		{
			_traceBack(particles[i], trace, pset);
		}
	}
	return trace;
}

void
_splitNow(vector<MovingParticle*>& points, //vectorized particles
		int beg, int end, //indices specifying the subsequence of points
		vector<int>& match, //if match[i]==j, then points[i] and points[j] overlaps
		vector<vector<MovingParticle*>>& polygons //add a region as being found.
		)
{
	if (beg >= end) return;
	vector<MovingParticle*> poly;
	for (int i = beg; i <= end;)
	{
		poly.push_back(points[i]);
		if (match[i] > i)
		{
			_splitNow(points, i + 1, match[i], match, polygons);
			i = match[i] + 1;
		}
		else
		{
			i++;
		}
	}
	if (poly.size() > 2)
	{
		polygons.push_back(poly);
	}
}

/*
From a vectorized particles, extract subsets that form closed regions.
*/
vector<vector<MovingParticle*>>
MovingParticle::closedRegions(vector<MovingParticle*>& points)
{
	vector<vector<MovingParticle*>> polygons;
	vector<int> match(points.size(), -1);
	float eps = 1.0e-5;
	for (int i = 0; i < points.size(); ++i)
	{
		int d = points.size();
		for (int j = 0; j < points.size(); ++j)
		{
			if (i == j) continue;
			if (Distance(points[i]->p0, points[j]->p0) < eps)
			{
				int k = Abs(i - j);
				if (k < d)
				{
					match[i] = j;
				}
			}
		}
	}

	_splitNow(points, 0, points.size() - 1, match, polygons);

	return polygons;
}

/*
From a vectorized particles, extract subsets that form closed regions IN TERMS OF INITIAL PARTICLES.
This is a faster version of closedRegion and does not work when we have to consider OFFSET POLYGONS.
*/
vector<vector<MovingParticle*>>
MovingParticle::closedRegions2(vector<MovingParticle*>& points)
{
	vector<vector<MovingParticle*>> polygons;
	set<StationaryParticle*> pset;
	map<StationaryParticle*, int> pmap;
	vector<MovingParticle*> vstack; //vector being used as a stack
	int first = 0, last = 0;
	for (int i = 0; i < points.size(); ++i)
	{
		StationaryParticle* sp = points[i]->init_particle;
		if (pset.find(sp) != pset.end())
		{
			vector<MovingParticle*> poly;
			int k = pmap[sp];
			if (k >= vstack.size())
			{
				printf("%d %d\n", i, k);
				for (int m = 0; m<points.size(); ++m)
				{
					printf("%d %f %f\n", points[m]->id, points[m]->p.m_X, points[m]->p.m_Y);
				}
				mexErrMsgTxt("Failed tracing in closedRegions2.\n");
			}
			poly.insert(poly.begin(), vstack.begin() + k, vstack.end());
			if (poly.size()>2)
			{
				polygons.push_back(poly);
			}
			vstack.erase(vstack.begin() + k, vstack.end());
		}
		else
		{
			pset.insert(sp);
			if (i == points.size() - 1) //wrap around the loop to end it
			{
				StationaryParticle* sp0 = points[0]->init_particle;
				vector<MovingParticle*> poly = vstack;
				if (poly.size() > 2)
				{
					poly.push_back(points[i]);
					polygons.push_back(poly);
				}
				vstack.clear();
			}
		}
		pmap[sp] = vstack.size(); //mark the last encounter of this.
		vstack.push_back(points[i]);
	}
	return polygons;
}
