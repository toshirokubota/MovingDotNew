#include <ParticleSimulatorGreedy.h>
#include <szmexutilitytemplate.h>
#include <IntersectionConvexPolygons.h>

/*
Find a precise location where the edge should be split by a split event.
*/
CParticleF
makeSplitParticle(EventStruct ev)
{
	CParticleF c = ev.p->move(ev.t);
	CParticleF q = ev.q->move(ev.t);
	CParticleF r = ev.r->move(ev.t);
	pair<float, float> param = _IntersectConvexPolygon::intersect(ev.p->getP(), c, q, r);
	float t = param.second;
	CParticleF q0 = ev.q->getP();
	CParticleF r0 = ev.r->getP();
	CParticleF f((1 - t)*q0.m_X + t*r0.m_X, (1 - t)*q0.m_Y + t*r0.m_Y);
	return f;
}

float
_fitnessMeasure(vector<CParticleF>& vp, float scale)
{
	float area = polygonArea(vp);
	float maxlen = 0;
	float sum = 0, sum2 = 0;
	for (int i = 0; i < vp.size(); ++i)
	{
		CParticleF p = vp[i];
		CParticleF q = vp[(i + 1) % vp.size()];
		//CParticleF r = vp[(i - 1 + vp.size()) % vp.size()];
		float d = Distance(p, q);
		maxlen = Max(d, maxlen);
		sum += d;
		sum2 += d * d;
	}
	//return sqrt(area) / maxlen;
	//return sqrt(area) / (scale + sum);
	//return area / (sum2);
	return sqrt(area) / scale;
}

/*
from a vector of particles (VP), find a polygon with non-zero area that contains the particle, P.
*/
vector<MovingParticle*>
extractSimplePath(MovingParticle* p0, MovingParticle* pend)
{
	struct _Linker {
		MovingParticle* p;
		_Linker* prev;
		_Linker(MovingParticle* p)
		{
			this->p = p;
			prev = NULL;
		}
	};

	int id = p0->getId(); //used as a unique number to detect looping

	vector<_Linker*> vl;
	_Linker* pl = new _Linker(p0);
	p0->getInitParticle()->link = (void*)pl;
	p0->getInitParticle()->info = id;
	vl.push_back(pl);
	MovingParticle* p = p0->getNext();
	while (p != pend)
	{
		if (p == p0) //could not reach pend. Cannot continue
		{
			mexErrMsgTxt("Failed to find a path from p0 to pend in ParticleSimulatorGreedy::extractClosedPolygon().");
		}
		if (p == NULL || pl == NULL) //could not reach pend. Cannot continue
		{
			mexErrMsgTxt("Null pointer in ParticleSimulatorGreedy::extractClosedPolygon().");
		}
		int id0 = p->getInitParticle()->info;
		_Linker* pl2 = new _Linker(p);
		if (id0 == id) //a loop is found. Break it.
		{
			_Linker* pl3 = (_Linker*)p->getInitParticle()->link;
			pl2->prev = pl3->prev;
		}
		else
		{
			pl2->prev = pl;
		}

		vl.push_back(pl2);
		pl = pl2;
		p->getInitParticle()->link = (void*)pl;
		p->getInitParticle()->info = id;
		p = p->getNext();
	}

	//trace back
	vector<MovingParticle*> area;
	while (pl != NULL)
	{
		area.insert(area.begin(), pl->p);
		pl = pl->prev;
	}
	//clean up
	for (int i = 0; i < vl.size(); ++i)
	{
		vl[i]->p->getInitParticle()->info = -1;
		vl[i]->p->getInitParticle()->link = NULL;

		delete vl[i];
		vl[i] = NULL;
	}

	return area;
}

float
ParticleSimulatorGreedy::computeFitness(vector<MovingParticle*>& area, float scale)
{
	vector<CParticleF> vp0;
	set<MovingParticle*> pset;
	for (int i = 0; i < area.size(); ++i)
	{
		vp0.push_back(area[i]->getP());
		pset.insert(area[i]);
	}
	float fit = _fitnessMeasure(vp0, scale);
	float coverage = 1.0f;
	for (int i = 0; i < chosen.size(); ++i)
	{
		set<MovingParticle*> iset;
		set<MovingParticle*> uset;
		set_intersection(pset.begin(), pset.end(),
			chosen[i]->getParticleSet().begin(), chosen[i]->getParticleSet().end(),
			std::inserter(iset, iset.begin()));
		set_union(pset.begin(), pset.end(),
			chosen[i]->getParticleSet().begin(), chosen[i]->getParticleSet().end(),
			std::inserter(uset, uset.begin()));
		coverage = Min(coverage, 1 - (float)iset.size() / (float)uset.size());
	}
	return fit * coverage;
}

pair<vector<MovingParticle*>,vector<MovingParticle*>>
split(EventStruct ev)
{
	pair<vector<MovingParticle*>, vector<MovingParticle*>> result;
	if (ev.type == SplitEvent || ev.type == CollisionEvent)
	{
		MovingParticle* pe = (MovingParticle*)ev.p;
		MovingParticle* qe = (MovingParticle*)ev.q;
		MovingParticle* re = (MovingParticle*)ev.r;
		if (pe != NULL && qe != NULL && re != NULL)
		{
			CParticleF f = makeSplitParticle(ev);
			{
				vector<MovingParticle*> area = extractSimplePath(pe, qe);
				if (area.empty() == false)
				{
					area.push_back(qe);
					result.first = area;
				}
			}
			{
				vector<MovingParticle*> area = extractSimplePath(pe, re);
				if (area.empty() == false)
				{
					area.push_back(re);
					result.second = area;
				}
			}
		}
	}
	return result;
}


bool
ParticleSimulatorGreedy::Simulate(float endtime, float delta, bool bdebug)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	PolygonFactory& pfactory = PolygonFactory::getInstance();
	float snapTime = delta;
	bool bSuccess = true;
	{
		vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
		snapshots.insert(snapshots.end(), shots.begin(), shots.end());
	}
	float thres = delta;
	int iter = 0;
	map<MovingParticle*, float> fmap;
	while (iter < endtime)
	{
		iter++;
		vector<pair<float, MovingParticle*>> pairs;
		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (p->isUnstable()) continue;
			if (p->isReflex())
			{
				float fitness = -std::numeric_limits<float>::infinity();
				if (fmap.find(p) != fmap.end())
				{
					fitness = fmap[p];
				}
				Polygon* polygon = p->getPolygon();
				EventStruct ev = p->getEvent();
				if (polygon->contains((MovingParticle*)ev.q) && polygon->contains((MovingParticle*)ev.r))
				{
				}
				else
				{
					p->setEvent(p->findNextSplitEvent(p->getPolygon()->getParticles()));
					pair<vector<MovingParticle*>, vector<MovingParticle*>> areas = split(p->getEvent());
					fitness1 = computeFitness(areas.first);
					if (fitness != fitness)
					{
						fitness = -std::numeric_limits<float>::infinity();
					}
				}
				pairs.push_back(pair<float, MovingParticle*>(fitness, p));
				fmap[p] = fitness;
			}
		}
		if (pairs.empty()) break;

		sort(pairs.begin(), pairs.end());
		MovingParticle* p = (*(pairs.end() - 1)).second;
		float fvalue = (*(pairs.end() - 1)).first;
		if (fvalue <= thres || fvalue != fvalue)
			break;

		//apply the event
		pair<MovingParticle*,MovingParticle*> pnew = applyEventGreedy(p->getEvent());
		{
			vector<Snapshot> shots = Snapshot::TakeSnapshot((float)iter);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			snapTime += delta;
		}
		printf("%d: fitness = %f\n", iter, fvalue);
		p->getEvent().print();

		MovingParticle* ap[2] = { pnew.first, pnew.second };
		for (int i = 0; i < 2; ++i)
		{
			if (ap[i] != NULL)
			{
				/*vector<vector<MovingParticle*>> areas = MovingParticle::closedRegions2(MovingParticle::vectorize(ap[i]));
				for (int j = 0; j < areas.size(); ++j)
				{
					Snapshot shot((float)iter, 0.0f, areas[j]);
					closedRegions.push_back(shot);
				}*/
				vector<MovingParticle*> area = extractSimplePath(ap[i], ap[i]);
				Snapshot shot((float)iter, 0.0f, area);
				closedRegions.push_back(shot);
				chosen.push_back(pfactory.makePolygon(area, (float)iter));
			}
		}

		//chosen.push_back(p->getPolygon());
		doneEvents.push_back(p->getEvent());
	}
	return bSuccess;
}

pair<MovingParticle*,MovingParticle*>
ParticleSimulatorGreedy::applyEventGreedy(EventStruct ev)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	pair<MovingParticle*, MovingParticle*> particles(NULL, NULL);
	if (ev.type == SplitEvent || ev.type == CollisionEvent)
	{
		MovingParticle* p = (MovingParticle*)ev.p;
		MovingParticle* q = (MovingParticle*)ev.q;
		MovingParticle* r = (MovingParticle*)ev.r;
		MovingParticle* pn = p->getNext();
		MovingParticle* pr = p->getPrev();
		StationaryParticle* sp = sfactory.makeParticle(makeSplitParticle(ev));
		{ //from q to p
			MovingParticle* pnew = factory.makeParticle(p->getInitParticle(), Split, 0.0f);
			MovingParticle* qnew = factory.makeParticle(sp, Split, 0.0f);
			ParticleDirection pd(sp->getP(), p->getP(), true);
			MovingFront mf(pd);
			MovingParticle::setNeighbors(qnew, q, pnew, q->getRear(), mf);
			MovingParticle::setNeighbors(pnew, qnew, pn, mf, p->getFront());
			qnew->calculateVelocityR();
			pnew->calculateVelocityR();
			MovingParticle::updatePolygon(pnew);
			particles.first = pnew;
		}
		{ //from p to r
			MovingParticle* pnew = factory.makeParticle(p->getInitParticle(), Split, 0.0f);
			MovingParticle* qnew = factory.makeParticle(sp, Split, 0.0f);
			ParticleDirection pd(p->getP(), sp->getP(), true);
			MovingFront mf(pd);
			MovingParticle::setNeighbors(qnew, pnew, r, mf, r->getRear());
			MovingParticle::setNeighbors(pnew, pr, qnew, p->getRear(), mf);
			qnew->calculateVelocityR();
			pnew->calculateVelocityR();
			MovingParticle::updatePolygon(pnew);
			particles.second = pnew;
		}
		factory.inactivate(p);
	}
	return particles;
}