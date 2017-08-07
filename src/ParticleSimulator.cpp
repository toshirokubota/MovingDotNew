#include <ParticleSimulator.h>
#include <map>
#include <szMiscOperations.h>
#include <MiscGeometry.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <szmexutilitytemplate.h>
#include <Polygon.h>

bool
ParticleSimulator::Simulate(float endtime, float delta, bool bdebug)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	PolygonFactory& pfactory = PolygonFactory::getInstance();
	float snapTime = delta;
	bool bSuccess = true;
	for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
	{
		(*it)->updateEvent();
	}
	{
		vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
		snapshots.insert(snapshots.end(), shots.begin(), shots.end());
	}
	int iter = 0;
	while (time < endtime)
	{
		iter++;
		MovingParticle* p = MovingParticle::getNextEvent();
		if (p == NULL) break;
		if (p->getEvent().t > endtime) break;
		if (p->id == 19 && p->event.q->id == 22)
			iter += 0;

		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			(*it)->update(p->getEvent().t - time);
		}
		time = p->getEvent().t;
		pair<MovingParticle*, MovingParticle*> pnew = p->applyEvent();
		/*if (p->event.type == SplitEvent || p->event.type == CollisionEvent)
		{
			if (p->polygon->getId() == p->event.q->polygon->getId())
			{
				//split the same polygon into two
				dag.add(p->polygon, pnew.first->polygon);
				dag.add(p->polygon, pnew.second->polygon);
			}
			else {
				//merge two polygons into one
				dag.add(pnew.first->polygon, p->polygon);
				dag.add(pnew.first->polygon, p->event.q->polygon);
			}
		}*/
		p->getEvent().print();

		//MovingParticle::removeUnstable();
		if (time >= snapTime)
		{
			vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			snapTime += delta;
		}
		MovingParticle::quickFinish();
		if (bdebug && MovingParticle::sanityCheck() == false)
		{
			printf("Violation of sanity check found at %f.\n", time);
			vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			bSuccess = false;
			break;
		}

		//find closed regions
		vector<vector<MovingParticle*>> regions = MovingParticle::clusterParticles();
		for (int i = 0; i < regions.size(); ++i)
		{
			vector<MovingParticle*> tr = MovingParticle::traceBackPolygon(regions[i]);
			vector<vector<MovingParticle*>> areas = MovingParticle::closedRegions(tr);

			Snapshot shot0(time, time, regions[i]);
			for (int j = 0; j < areas.size(); ++j)
			{
				Snapshot shot(time, 0.0f, areas[j]);
				if (find(closedRegions.begin(), closedRegions.end(), shot) == closedRegions.end())
				{
					closedRegions.push_back(shot);
					//Polygon* poly = pfactory.makePolygon(areas[j], time);
					traces.push_back(shot0);
				}
			}
			
			/*if (find(traces.begin(), traces.end(), shot0) == traces.end())
			{
				polygons.push_back(Snapshot(time, time, regions[i]));
			}*/
		}

		doneEvents.push_back(p->getEvent());
		//for (set<MovingParticle*>::iterator it = factory->updateQueue.begin(); it != factory->updateQueue.end(); ++it)
		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			(*it)->updateEvent();
		}
		factory.updateQueue.clear();
	}
	for (int i = 0; i < factory.particles.size(); ++i)
	{
		//factory.particles[i]->printParentTree("\t");
	}
	return bSuccess;
}

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
			float ang = GetVisualAngle2(u->key->getX(), u->key->getY(), 
				w->key->getX(), w->key->getY(), 
				v->key->getX(), v->key->getY());
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
ParticleSimulator::initializePolygon(vector<Edge<StationaryParticle*>*>& edges)
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

	PolygonFactory& polyFactory = PolygonFactory::getInstance();
	Polygon* polygon = polyFactory.makePolygon(particles, 0.0f);
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->setPolygon(polygon);
	}
	return particles;
}

/*
Trace a forrest and create a polygon for each tree.
*/
bool
traceForrest(vector<Vertex<StationaryParticle*>*>& forrest, vector<Edge<StationaryParticle*>*>& edges)
{
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
				vector<MovingParticle*> points = ParticleSimulator::initializePolygon(path);
				bDone = false;
				break;
			}
		}
		if (bDone)
		{
			break;
		}
	}
	return true;
}

bool
ParticleSimulator::Prepare(vector<StationaryParticle*>& points, vector<pair<int, int>>& E, float delta0)
{
	GraphFactory<StationaryParticle*>& factory = GraphFactory<StationaryParticle*>::GetInstance();
	ParticleFactory& pfactory = ParticleFactory::getInstance();

	vector<Vertex<StationaryParticle*>*> vertices;
	for (int i = 0; i < points.size(); ++i)
	{
		vertices.push_back(factory.makeVertex(points[i]));
	}
	vector<Edge<StationaryParticle*>*> edges;
	for (int i = 0; i < E.size(); ++i)
	{
		pair<int, int> idx = E[i];
		Edge<StationaryParticle*>* edge = factory.makeEdge(vertices[idx.first], vertices[idx.second], 1.0f);
		Edge<StationaryParticle*>* edge2 = factory.makeEdge(vertices[idx.second], vertices[idx.first], 1.0f);
		vertices[idx.first]->Add(edge);
		vertices[idx.second]->Add(edge2);
		edges.push_back(edge);
		edges.push_back(edge2);
		edge->u->key->neighbors.insert(edge->v->key);
		edge->v->key->neighbors.insert(edge->u->key);
	}

	traceForrest(vertices, edges);

	time = 0.0f;
	for (set<MovingParticle*>::iterator it = pfactory.activeSet.begin(); it != pfactory.activeSet.end(); ++it)
	{
		(*it)->update(delta0);
	}
	time += delta0;
	
	factory.Clean();
	return true;
}

#include <mexFileIO.h>
mxArray*
ParticleSimulator::SaveParticles()
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	const int dims[] = { factory.particles.size(), ParticleDumpSize };
	vector<float> F(dims[0] * dims[1]);
	for (int i = 0; i < dims[0]; ++i)
	{
		MovingParticle* p = factory.particles[i];
		vector<float> v = p->dump2vector();
		for (int j = 0; j < dims[1]; ++j)
		{
			SetData2(F, i, j, dims[0], dims[1], v[j]);
		}
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}

mxArray*
ParticleSimulator::SaveDoneEvents()
{
	const int dims[] = { doneEvents.size(), 5 };
	vector<float> F(dims[0] * dims[1]);
	for (int i = 0; i < doneEvents.size(); ++i)
	{
		EventStruct ev = doneEvents[i];
		SetData2(F, i, 0, dims[0], dims[1], (float)ev.p->id);
		SetData2(F, i, 1, dims[0], dims[1], (float)ev.q->id);
		SetData2(F, i, 2, dims[0], dims[1], (float)ev.r->id);
		SetData2(F, i, 0, dims[0], dims[1], (float)ev.t);
		SetData2(F, i, 0, dims[0], dims[1], (float)ev.type);
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}



bool
ParticleSimulator::LoadParticles(vector<float>& F, const int* dims)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	this->time = 0.0f;
	map<int, MovingParticle*> id2particle;
	vector<int> vchild1;
	vector<int> vchild2;
	vector<int> vparent1;
	vector<int> vparent2;
	vector<int> vprev;
	vector<int> vnext;
	vector<int> veq;
	vector<int> ver;
	for (int i = 0; i < dims[0]; ++i)
	{
		int id = (int)GetData2(F, i, 0, dims[0], dims[1], -1.0f);
		float x0 = GetData2(F, i, 1, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float y0 = GetData2(F, i, 2, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float x = GetData2(F, i, 3, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float y = GetData2(F, i, 4, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float vx = GetData2(F, i, 5, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float vy = GetData2(F, i, 6, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		int pid = (int)GetData2(F, i, 7, dims[0], dims[1], -1.0f);
		int nid = (int)GetData2(F, i, 8, dims[0], dims[1], -1.0f);
		MovingParticleType type = int2ParticleType((int)GetData2(F, i, 9, dims[0], dims[1], -1.0f));
		float created = GetData2(F, i, 10, dims[0], dims[1], 0.0f);
		float time = GetData2(F, i, 11, dims[0], dims[1], 0.0f);
		float ref = GetData2(F, i, 12, dims[0], dims[1], 0.0f);
		bool bActive = GetData2(F, i, 13, dims[0], dims[1], 0.0f)>0.0f ? true : false;
		bool bInitialized = (bool)GetData2(F, i, 14, dims[0], dims[1], 0.0f)>0.0f ? true : false;
		bool bUnstable = (bool)GetData2(F, i, 15, dims[0], dims[1], 0.0f)>0.0f ? true : false;
		int parent1 = (int)GetData2(F, i, 16, dims[0], dims[1], -1.0f);
		int parent2 = (int)GetData2(F, i, 17, dims[0], dims[1], -1.0f);
		EventType etype = int2EventType((int)GetData2(F, i, 18, dims[0], dims[1], 0.0f));
		int eq = (int)GetData2(F, i, 19, dims[0], dims[1], -1.0f);
		int er = (int)GetData2(F, i, 20, dims[0], dims[1], -1.0f);
		float etime = GetData2(F, i, 21, dims[0], dims[1], 0.0f);
		float frontx = GetData2(F, i, 22, dims[0], dims[1], 0.0f);
		float fronty = GetData2(F, i, 23, dims[0], dims[1], 0.0f);
		float frontsp = GetData2(F, i, 24, dims[0], dims[1], 0.0f);
		float rearx = GetData2(F, i, 22, dims[0], dims[1], 0.0f);
		float reary = GetData2(F, i, 23, dims[0], dims[1], 0.0f);
		float rearsp = GetData2(F, i, 24, dims[0], dims[1], 0.0f);
		MovingParticle* p = factory.makeParticle(sfactory.makeParticle(CParticleF(x0, y0)), type, created);
		p->p.m_X = x;
		p->p.m_Y = y;
		p->v[0] = vx;
		p->v[1] = vy;
		p->time = time;
		p->bActive = bActive;
		//p->bInitialized = bInitialized;
		p->bUnstable = bUnstable;
		p->event.type = etype;
		p->event.t = etime;
		p->reflexive = ref;
		vparent1.push_back(parent1);
		vparent2.push_back(parent2);
		p->event.p = p;
		vprev.push_back(pid);
		vnext.push_back(nid);
		veq.push_back(eq);
		ver.push_back(er);
		id2particle[p->id] = p;
		if (p->bActive==false)
		{
			factory.inactivate(p);
		}
		if (p->time > this->time)
		{
			this->time = p->time;
		}
		p->front = MovingFront(ParticleDirection(frontx, fronty), frontsp);
		p->rear = MovingFront(ParticleDirection(rearx, reary), rearsp);
	}
	//now  set prev, next, etc.
	for (int i = 0; i < factory.particles.size(); ++i)
	{
		MovingParticle* p = factory.particles[i];
		p->prev = vprev[i] >= 0 ? id2particle[vprev[i]] : NULL;
		p->next = vnext[i] >= 0 ? id2particle[vnext[i]] : NULL;
		p->event.q = veq[i] >= 0 ? id2particle[veq[i]] : NULL;
		p->event.r = ver[i] >= 0 ? id2particle[ver[i]] : NULL;
		p->parents[0] = vparent1[i] >= 0 ? id2particle[vparent1[i]] : NULL;
		p->parents[1] = vparent2[i] >= 0 ? id2particle[vparent2[i]] : NULL;
		if (p->id >= MovingParticle::_id)
		{
			MovingParticle::_id = p->id + 1;
		}
	}

	return true;
}

mxArray*
ParticleSimulator::SaveConvexity()
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	const int dims[] = { factory.particles.size(), 3 };
	vector<float> F(dims[0] * dims[1]);
	for (int i = 0; i < dims[0]; ++i)
	{
		MovingParticle* p = factory.particles[i];
		CParticleF pr = p->project(p->created + 0.1);
		SetData2(F, i, 0, dims[0], dims[1], pr.m_X);
		SetData2(F, i, 1, dims[0], dims[1], pr.m_Y);
		SetData2(F, i, 2, dims[0], dims[1], p->reflexive<=0 ? 0.0f: 1.0f);
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}

