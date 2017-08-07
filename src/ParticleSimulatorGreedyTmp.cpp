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

	vector<MovingParticle*> area;
	int id = p0->getId(); //used as a unique number to detect looping

	vector<_Linker*> vl;
	_Linker* pl = new _Linker(p0);
	p0->getInitParticle()->link = (void*)pl;
	p0->getInitParticle()->info = id;
	vl.push_back(pl);
	MovingParticle* p = p0->getNext();
	while (p != pend)
	{
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
ParticleSimulatorGreedy::_getCoverage(vector<MovingParticle*>& area)
{
	set<MovingParticle*> pset;
	for (int i = 0; i < area.size(); ++i)
	{
		pset.insert(area[i]);
	}
	set<MovingParticle*> iset;
	set<MovingParticle*> uset;
	set_intersection(pset.begin(), pset.end(),
		covered.begin(), covered.end(),
		std::inserter(iset, iset.begin()));
	set_union(pset.begin(), pset.end(),
		covered.begin(), covered.end(),
		std::inserter(uset, uset.begin()));
	float coverage = 1 - (float)iset.size() / (float)uset.size();
	return coverage;
}

pair<float, float> 
computeFitnessFromTempPoly(TempPolygons& poly)
{
	pair<float, float> result;
	float scale = Distance(poly.split, poly.ev.p->getP()) * 2;
	//first polygon
	{
		vector<CParticleF> vp;
		for (int i = 0; i < poly.poly1.size(); ++i)
		{
			vp.push_back(poly.poly1[i]->getP0());
		}
		vp.push_back(poly.split);
		result.first = _fitnessMeasure(vp, scale);
	}
	//second polygon
	{
		vector<CParticleF> vp;
		for (int i = 0; i < poly.poly2.size(); ++i)
		{
			vp.push_back(poly.poly2[i]->getP0());
		}
		vp.push_back(poly.split);
		result.first = _fitnessMeasure(vp, scale);
	}
	return result;
}

/*FitnessStruct
ParticleSimulatorGreedy::computeFitness(EventStruct ev, bool left)
{
	FitnessStruct fs;
	MovingParticle* pe = (MovingParticle*)ev.p;
	MovingParticle* qe = (MovingParticle*)ev.q;
	MovingParticle* re = (MovingParticle*)ev.r;
	if (pe != NULL && qe != NULL && re != NULL)
	{
		CParticleF f = makeSplitParticle(ev);
		vector<MovingParticle*> area;
		if (left)
		{
			area = extractSimplePath(pe, qe);
			area.push_back(qe);
		}
		else
		{
			area = extractSimplePath(re, pe);
			area.insert(area.begin(), pe);
		}
		vector<CParticleF> vp0;
		for (int i = 0; i < area.size(); ++i)
		{
			vp0.push_back(area[i]->getP());
		}
		vp0.push_back(f);
		float coverage = _getCoverage(area);
		float fit = _fitnessMeasure(vp0, Distance(f, pe->getP()) * 2);
		float value = coverage * fit;
		if (value > fs.value)
		{
			fs.value = value;
			fs.fitness = fit;
			fs.coverage = coverage;
			fs.f = f;
			fs.area = area;
			fs.bleft = left;
		}
	}
	return fs;
}*/

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
	while (iter < endtime)
	{
		iter++;
		FitnessStruct bestFitness;
		MovingParticle* psel = NULL;
		bool bLeft = true;
		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (p->isUnstable()) continue;
			if (p->isReflex())
			{
				p->setEvent(p->findNextSplitEvent(p->getPolygon()->getParticles()));
				EventStruct ev = p->getEvent();
				TempPolygons tpoly = testEventGreedy(ev);
				pair<float, float> fit = computeFitnessFromTempPoly(tpoly);
				//FitnessStruct fitLeft = computeFitness(ev, true);
				//FitnessStruct fitRight = computeFitness(ev, false);
				if (fit.first > bestFitness.value)
				{
					bestFitness = fit.first;
					psel = p;
					bLeft = true;
				}
				if (fitRight.value > bestFitness.value)
				{
					bestFitness = fitRight;
					psel = p;
					bLeft = false;
				}*/
			}
		}
		MovingParticle* p = psel;
		float fvalue = bestFitness.value;
		if (fvalue <= thres || fvalue != fvalue)
			break;

		//apply the event
		pair<MovingParticle*,MovingParticle*> pnew = applyEventGreedy(p->getEvent());
		{
			vector<Snapshot> shots = Snapshot::TakeSnapshot((float)iter);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			snapTime += delta;
		}
		p->getEvent().print();

		{
			MovingParticle* ap;
			MovingParticle* ap2;
			if (bLeft)
			{
				ap = pnew.first;
				ap2 = pnew.second;
			}
			else
			{
				ap = pnew.second;
				ap2 = pnew.first;
			}
			vector<MovingParticle*> area = extractSimplePath(ap, ap);
			vector<MovingParticle*> area2 = extractSimplePath(ap2, ap2);
			Snapshot shot((float)iter, 0.0f, area);
			closedRegions.push_back(shot);
			chosen.push_back(pfactory.makePolygon(area, (float)iter));
			for (int j = 0; j < area.size(); ++j)
			{
				covered.insert(area[j]);
			}
			printf("%d: fitness = %f, %f, %f, %d, %d %s\n", 
				iter, bestFitness.value, bestFitness.fitness, bestFitness.coverage, bestFitness.area.size(), area.size(),
				bestFitness.bleft ? "Left": "Right");
			int dc = bestFitness.area.size() - area.size();
			if (Abs(dc)>=2)
			{
				printf("New particles:\n");
				pnew.first->print();
				pnew.first->getPrev()->print();
				pnew.second->print();
				pnew.second->getNext()->print();
				printf("Best fitness:\n");
				for (int j = 0; j < bestFitness.area.size(); ++j)
				{
					printf("%d %f %f\n", bestFitness.area[j]->getId(), bestFitness.area[j]->getP0().m_X, bestFitness.area[j]->getP0().m_Y);
				}
				printf("Area 1:\n");
				for (int j = 0; j < area.size(); ++j)
				{
					printf("%d %f %f\n", area[j]->getId(), area[j]->getP0().m_X, area[j]->getP0().m_Y);
				}
				printf("Area 2:\n");
				for (int j = 0; j < area2.size(); ++j)
				{
					printf("%d %f %f\n", area2[j]->getId(), area2[j]->getP0().m_X, area2[j]->getP0().m_Y);
				}
				vector<MovingParticle*> area = extractSimplePath(ap, ap);
				vector<MovingParticle*> area2 = extractSimplePath(ap2, ap2);
				mexErrMsgTxt("Premature exist.");
			}
		}

		//chosen.push_back(p->getPolygon());
		doneEvents.push_back(p->getEvent());
	}
	return bSuccess;
}

pair<MovingParticle*,MovingParticle*>
ParticleSimulatorGreedy::applyEventGreedy(TempPolygons& t)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	pair<MovingParticle*, MovingParticle*> particles(NULL, NULL);
	EventStruct ev = t.ev;
	if (ev.type == SplitEvent || ev.type == CollisionEvent)
	{
		MovingParticle* p = (MovingParticle*)ev.p;
		MovingParticle* q = (MovingParticle*)ev.q;
		MovingParticle* r = (MovingParticle*)ev.r;
		MovingParticle* pr = p->getPrev();
		MovingParticle* pn = p->getNext();
		StationaryParticle* sp = sfactory.makeParticle(t.split);
		StationaryParticle* pi = p->getInitParticle();
		{ //Polygon 1
			MovingParticle* pnew = factory.makeParticle(pi, Split, 0.0f);
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

TempPolygons
ParticleSimulatorGreedy::testEventGreedy(EventStruct ev)
{
	TempPolygons tpoly;
	if (ev.type == SplitEvent || ev.type == CollisionEvent)
	{
		MovingParticle* p = (MovingParticle*)ev.p;
		MovingParticle* q = (MovingParticle*)ev.q;
		MovingParticle* r = (MovingParticle*)ev.r;
		MovingParticle* pn = p->getNext();
		MovingParticle* pr = p->getPrev();
		tpoly.split = makeSplitParticle(ev);

		{ //from p to q
			tpoly.poly1 = MovingParticle::extractPath(p, q);
		}
		{ //from p to r
			tpoly.poly2 = MovingParticle::extractPath(r, p);
			//move p to the beginning
			tpoly.poly2.insert(tpoly.poly2.begin(), p);
			tpoly.poly2.erase(tpoly.poly2.end() - 1);
		}
	}
	return tpoly;
}
