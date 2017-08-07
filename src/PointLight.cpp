#include <PointLight.h>
#include <szMiscOperations.h>
using namespace _PointLightNS;

float
_Distance(PointLight* p, PointLight* q)
{
	float a = 1.0f;
	return sqrt((p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y) + a * (p->frame - q->frame)*(p->frame - q->frame));
}

/*
A stationary candidate has all tp.p, tp.q, and tp.r are all same.
*/
bool
PointLight::isStationary(Candidate& c)
{
	return c.tp.p == c.tp.q && c.tp.q == c.tp.r;
}


/*
Measure how p is compatible to q.
This assumes causality, thus p has to be later frame than q to be compatible.
*/
float
PointLight::_Compatibility(Candidate& a, Candidate& b)
{
	float sgm = 10.0f;
	if (a.tp.q == b.tp.q) return Threshold; //this happens for points in the first or last frames as well as isolated points.
	float val = 0;
	if (a.tp.p == b.tp.q && a.tp.q == b.tp.r)
	{
		CParticleF a0 = a.tp.evaluate(a.tp.q->frame);
		CParticleF a1 = a.tp.evaluate(a.tp.p->frame);
		CParticleF b0 = b.tp.evaluate(a.tp.q->frame);
		CParticleF b1 = b.tp.evaluate(a.tp.p->frame);
		float d1 = Distance(a0, b0);
		float d2 = Distance(a1, b1);
		CParticleF u0 = a.tp.velocity(a.tp.q->frame);
		CParticleF u1 = a.tp.velocity(a.tp.p->frame);
		CParticleF v0 = b.tp.velocity(a.tp.q->frame);
		CParticleF v1 = b.tp.velocity(a.tp.p->frame);
		float d3 = Distance(u0, v0);
		float d4 = Distance(u1, v1);
		val = exp(-(d1*d1 + d2*d2 + d3*d3 + d4*d4) / (2 * sgm*sgm));
	}
	else if (a.tp.r == b.tp.q && a.tp.q == b.tp.p)
	{
		CParticleF a0 = a.tp.evaluate(a.tp.q->frame);
		CParticleF a1 = a.tp.evaluate(a.tp.r->frame);
		CParticleF b0 = b.tp.evaluate(a.tp.q->frame);
		CParticleF b1 = b.tp.evaluate(a.tp.r->frame);
		float d1 = Distance(a0, b0);
		float d2 = Distance(a1, b1);
		CParticleF u0 = a.tp.velocity(a.tp.q->frame);
		CParticleF u1 = a.tp.velocity(a.tp.r->frame);
		CParticleF v0 = b.tp.velocity(a.tp.q->frame);
		CParticleF v1 = b.tp.velocity(a.tp.r->frame);
		float d3 = Distance(u0, v0);
		float d4 = Distance(u1, v1);
		val = exp(-(d1*d1 + d2*d2 + d3*d3 + d4*d4) / (2 * sgm*sgm));
	}
	else
	{
		val = 0.0f;
	}
	return val;
}

float
PointLight::_CompatibilityStationary(Candidate& a, Candidate& b)
{
	float sgm = StationarySigma;
	float d = _Distance(a.tp.q, b.tp.q);
	float val = exp(-d*d / (2 * sgm*sgm));
	return val;
}

void
PointLight::print(char* tab, char* newl)
{
	if (tab != NULL) printf("%s", tab);
	printf("%d(%3.3f,%3.3f,%f) [", id, x, y, frame);
	if (newl != NULL) printf("%s", newl);
}

void
PointLight::initializeProb()
{
	int n = candidates.size();
	for (int i = 0; i < candidates.size(); ++i)
	{
		candidates[i].prob = 1.0 / n;
	}
}

PointLight::PointLightTriple 
PointLight::winner()
{ 
	float maxP = 0;
	float totalP = 0;
	int idx = -1;
	for (int i = 0; i < candidates.size(); ++i)
	{
		if (candidates[i].prob > maxP)
		{
			maxP = candidates[i].prob;
			idx = i;
		}
		totalP += candidates[i].prob;
	}
	float eta = 1.0 - totalP;
	if (eta > maxP || idx < 0)
	{
		return PointLightTriple(NULL, this, NULL);
	}
	else
	{
		return candidates[idx].tp;
	}
}

/*
ETA is a common value added to each candidate. This is to slow down the committment so that the overall configuration does not
quickly converge to a local optimum.
*/
void
PointLight::updateFitness()
{
	if (bDirty == false) return;
	float rate = 0.;

	for (int i = 0; i < candidates.size(); ++i)
	{
		if (isStationary(candidates[i]))
		{
			candidates[i].comp = Threshold;
		}
		else
		{
			float s1 = Threshold;
			for (int j = 0; j < candidates[i].tp.p->candidates.size(); ++j)
			{
				float fit = _Compatibility(candidates[i], candidates[i].tp.p->candidates[j]) + candidates[i].tp.p->candidates[j].comp * rate;
				s1 += fit * candidates[i].tp.p->candidates[j].prob;
				//s = Max(s, _Compatibility(candidates[i], candidates[i].tp.p->candidates[j]) * candidates[i].tp.p->candidates[j].prob);
				if (id == -1)
				{
					printf(">(%d %d %d),(%d %d %d):: %f, %f, %f, %f\n", 
						candidates[i].tp.p->id, candidates[i].tp.q->id, candidates[i].tp.r->id,
						candidates[i].tp.p->candidates[j].tp.p->id, candidates[i].tp.p->candidates[j].tp.q->id, candidates[i].tp.p->candidates[j].tp.r->id,
						_Compatibility(candidates[i], candidates[i].tp.p->candidates[j]), candidates[i].tp.p->candidates[j].prob, fit, s1);
				}
			}
			float s2 = Threshold;
			for (int j = 0; j < candidates[i].tp.r->candidates.size(); ++j)
			{
				float fit = _Compatibility(candidates[i], candidates[i].tp.r->candidates[j]) + candidates[i].tp.r->candidates[j].comp * rate;
				s2 += fit * candidates[i].tp.r->candidates[j].prob;
				if (id == -1)
				{
					printf(">(%d %d %d),(%d %d %d):: %f, %f, %f, %f\n", 
						candidates[i].tp.p->id, candidates[i].tp.q->id, candidates[i].tp.r->id,
						candidates[i].tp.r->candidates[j].tp.p->id, candidates[i].tp.r->candidates[j].tp.q->id, candidates[i].tp.r->candidates[j].tp.r->id,
						_Compatibility(candidates[i], candidates[i].tp.r->candidates[j]), candidates[i].tp.r->candidates[j].prob, fit, s2);
				}
				//s = Max(s, _Compatibility(candidates[i], candidates[i].tp.r->candidates[j]) * candidates[i].tp.r->candidates[j].prob);
			}
			candidates[i].comp = sqrt(s1 * s2);
			if (id == -1)
			{
				//if (candidates[i].comp > 0.1)
				{
					printf("F %d) %d (%d %d %d) => %f (%f)\n", id, i, candidates[i].tp.p->id, candidates[i].tp.q->id, candidates[i].tp.r->id,
					candidates[i].comp, candidates[i].prob);
				}
			}
		}
	}
}

/*
ETA is a common value added to each candidate. This is to slow down the committment so that the overall configuration does not
quickly converge to a local optimum.
*/
void
PointLight::updateProb()
{
	if (bDirty == false) return;

	float totalP = 0;
	float total = 0;
	for (int i = 0; i < candidates.size(); ++i)
	{
		total += candidates[i].comp * candidates[i].prob;
		totalP += candidates[i].prob;
	}
	//float eta = Threshold; // totalP >= 1.0 ? 0.0f : (1.0 - totalP) * Threshold;
	float eta = 0.0f; // Threshold;
	for (int i = 0; i < candidates.size(); ++i)
	{
		candidates[i].prob = candidates[i].comp * candidates[i].prob / (total + eta);
		/*if (id == 39 || id == 20 || id == 30 || id == 29) // || id==41 || id==101)
		{
			if (candidates[i].prob > 0.1)
			{
				printf("P %d) %d (%d %d %d) => %f (%f)\n", id, i, candidates[i].tp.p->id, candidates[i].tp.q->id, candidates[i].tp.r->id, 
					candidates[i].comp, candidates[i].prob);
			}
		}*/
	}
}

void
PointLight::updateState()
{
	float maxp = 0;
	int idx = -1;
	for (int i = 0; i < candidates.size(); ++i)
	{
		if (candidates[i].prob > maxp)
		{
			maxp = candidates[i].prob;
			idx = i;
		}
	}
	if (maxp <= Threshold) state = Isolated;
	else {
		PointLightTriple tp = candidates[idx].tp;
		if (tp.p == tp.q && tp.r == tp.q) state = Isolated;
		else if (tp.p == tp.q) state = Appeared;
		else if (tp.r == tp.q) state = Seeking;
		else state = Connected;
	}
	bDirty = false;
	if (id == 403 || id == 404)
	{
		if (idx >= 0)
		{
			printf("%d %d, winner=(%d,%d,%d)\n", id, state, candidates[idx].tp.p->id, candidates[idx].tp.q->id, candidates[idx].tp.r->id);
		}
		else
		{
			printf("%d %d, winner=NA\n", id, state);
		}
	}
}

#include <set>
void
PointLight::updateCandidates() {
	if (candidates0.empty() == false)
	{
		bDirty = true;
		for (int i = 0; i < candidates0.size(); ++i)
		{
			if (find(candidates.begin(), candidates.end(), candidates0[i]) == candidates.end())
			{
				candidates.push_back(candidates0[i]);
			}
		}
		/*for (int i = 0; i < candidates.size(); ++i)
		{
			candidates[i].tp.p->bDirty = true;
			candidates[i].tp.r->bDirty = true;
		}*/
		candidates0.clear();
		initializeProb();
	}
}

/*
A convenient method to calculate an intermediate 6-vector for EQW update.
The vector will need to be multiplied with G vector.
*/
void
_CalculateEQWPairTerm(PointLight::Candidate& a, PointLight::Candidate& b, float res[6], float H[3][3])
{
	PointLight*  p = b.tp.q;
	if (b.tp.q->frame > a.tp.q->frame) {
		for (int i = 0; i < 3; ++i)
		{
			res[i] = H[i][0] * b.tp.v[0] + H[i][1] * b.tp.v[1] + H[i][2] * b.tp.v[2];
		}
		for (int i = 0; i < 3; ++i)
		{
			res[i+3] = H[i][0] * b.tp.v[3] + H[i][1] * b.tp.v[4] + H[i][2] * b.tp.v[5];
		}
	}
	else
	{
		for (int i = 0; i < 3; ++i)
		{
			res[i] = H[0][i] * b.tp.v[0] + H[1][i] * b.tp.v[1] + H[2][i] * b.tp.v[2];
		}
		for (int i = 0; i < 3; ++i)
		{
			res[i + 3] = H[0][i] * b.tp.v[3] + H[1][i] * b.tp.v[4] + H[2][i] * b.tp.v[5];
		}
	}
}

void
PointLight::updateParams()
{
	float r = 1.0f, t = 1.0f;
	float a = 4 * r + t*t, b = 8 * r + t*t, c = 2 * r + t*t;
	float G[3][3], H[3][3];
	G[0][0] = .5*a / b; G[0][1] = -.5 / b; G[0][2] = 0;
	G[1][0] = -.5 / b; G[1][1] = 1. / (t*t*b); G[1][2] = 0;
	G[2][0] = 0; G[2][1] = 0; G[2][2] = .5 / c;
	H[0][0] = -2.; H[0][1] = -t*t; H[0][2] = t;
	H[1][0] = -t*t; H[1][1] = 0; H[1][2] = -2.*r*t;
	H[2][0] = -t; H[2][1] = 2.*r*t; H[2][2] = -2.*r;

	for (int i = 0; i < candidates.size(); ++i)
	{
		Candidate a = candidates[i];
		PointLight* p = a.tp.p;
		PointLight* q = a.tp.q;
		PointLight* r = a.tp.r;
		if (p == q || p == r) continue;
		float res[6] = { 0, 0, 0, 0, 0, 0 };
		for (int j = 0; j < p->candidates.size(); ++j)
		{
			Candidate b = p->candidates[j];
			if (a.tp.p == b.tp.q && a.tp.q == b.tp.r)
			{
				float res0[6];
				_CalculateEQWPairTerm(a, b, res0, H);
				for (int k = 0; k < 6; ++k)
				{
					res[k] += b.prob * res0[k];
				}
			}
		}
		for (int j = 0; j < r->candidates.size(); ++j)
		{
			Candidate b = r->candidates[j];
			if (a.tp.r == b.tp.q && a.tp.q == b.tp.p)
			{
				float res0[6];
				_CalculateEQWPairTerm(a, b, res0, H);
				for (int k = 0; k < 6; ++k)
				{
					res[k] += b.prob * res0[k];
				}
			}
		}
		for (int k = 0; k < 3; ++k)
		{
			candidates[i].tp.v[k] = G[k][0] * res[0] + G[k][1] * res[1] + G[k][2] * res[2];
			candidates[i].tp.v[k+3] = G[k][0] * res[3] + G[k][1] * res[4] + G[k][2] * res[5];
		}
	}
}