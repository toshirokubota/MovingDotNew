#include <OffsetPolygonDAGBuilder.h>
#include <map>
#include <szMiscOperations.h>
#include <MiscGeometry.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <szmexutilitytemplate.h>
#include <IntersectionConvexPolygons.h>

vector<Edge<StationaryParticle*>*>
_trace(Edge<StationaryParticle*>* edge)
{
	Vertex<StationaryParticle*>*first = edge->u;
	vector<Edge<StationaryParticle*>*> path;

	edge->type = Back; //mark it Back once visited.
	path.push_back(edge);
	while (true)
	{
		//choose the first one in the clockwise order
		Edge<StationaryParticle*>* next = NULL;
		float minAngle = 3 * PI;
		//Walk in clock-wise  order.
		Vertex<StationaryParticle*>* u = edge->u;
		Vertex<StationaryParticle*>* v = edge->v;
		for (int i = 0; i < v->aList.size(); ++i)
		{
			Vertex<StationaryParticle*>* w = v->aList[i]->v;
			float ang = GetVisualAngle2(u->key->getP().m_X, u->key->getP().m_Y, w->key->getP().m_X, w->key->getP().m_Y, v->key->getP().m_X, v->key->getP().m_Y);
			if (ang <= 0) ang += 2 * PI;
			if (ang < minAngle)
			{
				minAngle = ang;
				next = v->aList[i];
			}
		}
		if (next->type == Back) break; //Done.

		edge = next;
		edge->type = Back; //mark it Back once visited.
		path.push_back(edge);
	}

	return path;
}

vector<MovingParticle*>
initializePolygon(vector<Edge<StationaryParticle*>*>& edges)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<MovingParticle*> particles;
	for (int i = 0; i < edges.size(); ++i)
	{
		particles.push_back(factory.makeParticle(edges[i]->u->key, Initial, 0.0f));
		if (edges[i]->u->aList.size() <= 1) //leaf node.
		{
			particles.push_back(factory.makeParticle(edges[i]->u->key, Initial, 0.0f));
		}
	}
	//set the neighborhood
	for (int i = 0; i < particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		MovingParticle* q = particles[i == particles.size() - 1 ? 0 : i + 1];
		MovingParticle* r = particles[i == 0 ? particles.size() - 1 : i - 1];
		MovingParticle::setNeighbors(p, r, q);
	}
	//calculate velocity 
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->initializeVelocity();
	}
	return particles;
}

/*
Trace a forrest and create a polygon for each tree.
*/
vector<vector<MovingParticle*>>
traceForrest(vector<Edge<StationaryParticle*>*>& edges)
{
	vector<vector<MovingParticle*>> forrest;
	for (int i = 0; i < edges.size(); ++i)
	{
		edges[i]->type = Tr;
	}

	while (true)
	{

		bool bDone = true;
		for (int j = 0; j < edges.size(); ++j)
		{
			if (edges[j]->type == Tr)
			{
				vector<Edge<StationaryParticle*>*> path = _trace(edges[j]);
				vector<MovingParticle*> points = initializePolygon(path);
				forrest.push_back(points);
				bDone = false;
				break;
			}
		}
		if (bDone)
		{
			break;
		}
	}
	return forrest;
}

bool
OffsetPolygonDAGBuilder::Polygonify()
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	factory.clean();

	forrest = traceForrest(dag.edges);

	float delta = 0.01;
	for (int i = 0; i < forrest.size(); ++i)
	{
		for (int j = 0; j < forrest[i].size(); ++j)
		{
			MovingParticle* p = forrest[i][j];
			p->update(delta);
			StationaryParticle* sp = p->getInitParticle();
			if (p->isLeaf() == false) aMap[sp] = false;
		}
	}
	for (int i = 0; i < forrest.size(); ++i)
	{
		for (int j = 0; j < forrest[i].size(); ++j)
		{
			MovingParticle* p = forrest[i][j];
			p->updateEvent();
		}
	}
	return true;
}

bool
OffsetPolygonDAGBuilder::Expand(float endtime)
{
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	ParticleFactory& pfactory = ParticleFactory::getInstance();
	float time = 0;
	float cutoff = 0.5;
	while (time <= endtime)
	{
		MovingParticle* p = MovingParticle::getNextEvent();
		if (p == NULL) break;

		if (aMap[p->getInitParticle()] //still active
			&& (p->event.type == CollisionEvent || p->event.type == SplitEvent)) //->isReflex()) 
		{
			MovingParticle* q = (MovingParticle*)p->event.q;
			MovingParticle* r = (MovingParticle*)p->event.r;
			CParticleF pp = p->move(p->event.t);
			CParticleF qp = q->move(p->event.t);
			CParticleF rp = r->move(p->event.t);
			pair<float,float> param = _IntersectConvexPolygon::intersect(p->getP0(), pp, qp, rp);
			if (param.second < cutoff)
			{
				dag.connect(p->getInitParticle(), q->getInitParticle(), time);
			}
			else 
			{
				dag.connect(p->getInitParticle(), r->getInitParticle(), time);
			}
			aMap[p->getInitParticle()] = false;
			aMap[q->getInitParticle()] = false;
		}
		time = p->event.t;
		pfactory.inactivate(p);
		
	}
	return true;
}
