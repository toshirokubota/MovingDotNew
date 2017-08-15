/*
This routine implements grouping of moving points.
*/
#include <mex.h>

#include <iostream>
#include <fstream>
#include <set>
#include <map>
using namespace std;
#include <stdlib.h>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>
#include <szDistanceTransform.h>
#include <szDistanceTransformNonIsotropic.h>
#include <MovingDot.h>
#include <Graph.h>
#include <CoreParticle.h>
#include <NeighborhoodFactory.h>
#include <DisjointSet.h>
#include <szMyNeighborOp.h>
#include <GraphFactory.h>
#include <Kruskal.h>

bool _EightNeighbor = false;

float
dotProduct(float x1, float y1, float z1, float t1, float x2, float y2, float z2, float t2)
{
	return x1*x2 + y1*y2 + z1*z2 + t1*t2;
}

/*
When two particles (p and q) meet originated from different strong core sduring the growth, and a graph edge will be formed.
This routine provides the weight of the edge.
*/
float weight(CoreParticle* p, CoreParticle* q)
{
	//return sqrt(p->src->dval * q->src->dval) / sqrt(p->dval * q->dval);
	//float len = length(p->src->x - q->src->x, p->src->y - q->src->y, p->src->z - q->src->z, p->src->t - q->src->t);
	//return (len + sqrt(p->src->dval * q->src->dval) + 1.0) / sqrt(p->dval * q->dval);
	float len = length(p->src->x - q->src->x, p->src->y - q->src->y, p->src->z - q->src->z, p->src->t - q->src->t);
	return len*len;
}

bool
surfaceParticle(CoreParticle* p, int ndim, const int* dims)
{
	if (ndim == 2) //08/08/2017
	{
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1)
		{
			return false;  //consider image to be continuing beyond the boundary.
		}
	}
	else if (ndim == 3) //08/10/2017
	{
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1 || p->z == 0 || p->z == dims[2]-1)
		{
			return false;  //consider image to be continuing beyond the boundary.
		}
	}
	else if (ndim == 4) //08/10/2017
	{
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1 || p->z == 0 || p->z == dims[2] - 1 || p->t == 0 || p->t == dims[3] - 1)
		{
			return false;  //consider image to be continuing beyond the boundary.
		}
	}
	int ns = NeighborhoodFactory::getInstance(ndim).neighbor4.size();
	if (_EightNeighbor)
	{
		ns = NeighborhoodFactory::getInstance(ndim).neighbor8.size();
	}
	return p->neighbors.size() < ns;
}

bool
medialParticle(CoreParticle* p, int ndim)
{
	return p->ascendents.size() >= 2 * (ndim - 1) || p->descendents.size() == 0;
}

vector<vector<CoreParticle*>>
clusterParticles(vector<CoreParticle*>& particles);

/*
Ascendents of p are not linearly separable.
*/
bool
strongMedialParticle(CoreParticle* p, int ndim, const int* dims)
{
	if (p->core >= 0) return p->core;
	//if (ndim == 2) //08/08/2017
	{
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1)
		{
			p->core = 0;
			return false;  //consider image to be continuing beyond the boundary.
		}
	}
	/*else if (ndim == 3) //08/10/2017
	{
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1 || p->z == 0 || p->z == dims[2] - 1)
		{
			p->core = 0;
			return false;  //consider image to be continuing beyond the boundary.
		}
	}
	else if (ndim == 4) //08/10/2017
	{
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1 || p->z == 0 || p->z == dims[2] - 1 || p->t == 0 || p->t == dims[3] - 1)
		{
			p->core = 0;
			return false;  //consider image to be continuing beyond the boundary.
		}
	}*/
	if (p->descendents.size() == 0)
	{
		p->core = 2;
		//return true;
	}
	else
	{
		p->core = 0;
	}
	if (ndim == 2 || ndim == 3)
	{
		/*if (p->ascendents.size() >= 2 * (ndim - 1))
		{
			p->core = 1;
		}*/
		/*if (p->dval >= 7.0f)
		{
		bool blmax = true;
		for (int i = 0; i<p->neighbors8.size(); ++i)
		{
		if (p->neighbors8[i]->dval >= p->dval)
		{
		blmax = false;
		break;
		}
		}
		if (blmax) p->core = 1;
		}*/
	}
	return p->core > 0;
}

void
labelStrongCoreParticles(vector<CoreParticle*>& particles, int ndim, const int* dims)
{
	/*vector<CoreParticle*> surfaces;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (surfaceParticle(particles[i], ndim))
		{
			surfaces.push_back(particles[i]);
		}
	}
	vector<vector<CoreParticle*>> groups = clusterParticles(surfaces);
	map<CoreParticle*, int> imap;*/
	for (int i = 0; i < particles.size(); ++i)
	{
		strongMedialParticle(particles[i], ndim, dims);
	}
}

bool
inflectionParticle(CoreParticle* p, int ndim)
{
	//return surfaceParticle(p, ndim) && p->descendents.size() >= 2;  //ndim;
	return false; //disabled -- 08/08/2017
}//

CoreParticle*
coreParticleNdim(vector<int>& loc, int ndim)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	switch (ndim)
	{
	case 1:
		return factory.makeParticle(loc[0]);
	case 2:
		return factory.makeParticle(loc[0], loc[1]);
	case 3:
		return factory.makeParticle(loc[0], loc[1], loc[2]);
	case 4:
		return factory.makeParticle(loc[0], loc[1], loc[2], loc[3]);
	default:
		mexErrMsgTxt("particleNdim: unsupported number of dimensions. It has to be between 1 and 4.");
		return NULL;
	}
}

vector<int>
coreParticle2Index(CoreParticle*  p, int ndim)
{
	vector<int> idx(ndim);
	if (ndim == 1)
	{
		idx[0] = p->x;
	}
	else if (ndim == 2)
	{
		idx[0] = p->x;
		idx[1] = p->y;
	}
	else if (ndim == 3)
	{
		idx[0] = p->x;
		idx[1] = p->y;
		idx[2] = p->z;
	}
	else if (ndim == 4)
	{
		idx[0] = p->x;
		idx[1] = p->y;
		idx[2] = p->z;
		idx[3] = p->t;
	}
	else
	{
		mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
	}
	return idx;
}

template<class T>
bool
SetVoxel(vector<T>& A,
	const CoreParticle* p,
	T value,
	int ndim,
	const int* dims)
{
	if (ndim == 1)
	{
		return SetData(A, p->x, value);
	}
	else if (ndim == 2)
	{
		return SetData2(A, p->x, p->y, dims[0], dims[1], value);
	}
	else if (ndim == 3)
	{
		return SetData3(A, p->x, p->y, p->z, dims[0], dims[1], dims[2], value);
	}
	else if (ndim == 4)
	{
		return SetData4(A, p->x, p->y, p->z, p->t, dims[0], dims[1], dims[2], dims[3], value);
	}
	else
	{
		mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
		return false;
	}
}

template<class T>
T
GetVoxel(const vector<T>& A,
	const CoreParticle* p,
	T defaultValue,
	int ndim,
	const int* dims)
{
	if (ndim == 1)
	{
		return GetData(A, p->x, defaultValue);
	}
	else if (ndim == 2)
	{
		return GetData2(A, p->x, p->y, dims[0], dims[1], defaultValue);
	}
	else if (ndim == 3)
	{
		return GetData3(A, p->x, p->y, p->z, dims[0], dims[1], dims[2], defaultValue);
	}
	else if (ndim == 4)
	{
		return GetData4(A, p->x, p->y, p->z, p->t, dims[0], dims[1], dims[2], dims[3], defaultValue);
	}
	else
	{
		mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
		return defaultValue;
	}
}

vector<CoreParticle*>
generateParticleMap(vector<unsigned char>& L,
	int ndim,
	const int* dims)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<CoreParticle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	for (int i = 0; i < L.size(); ++i)
	{
		if (L[i])
		{
			vector<int> sub = Ind2Sub(i, ndim, dims);
			CoreParticle* p = coreParticleNdim(sub, ndim);
			mp[i] = p;
		}
	}
	return mp;
}

vector<CoreParticle*>
setupParticleNeighbors(vector<CoreParticle*>& mp,
	int ndim,
	const int* dims)
{
	vector<CoreParticle*> particles;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i] != NULL)
		{
			particles.push_back(mp[i]);
		}
	}
	NeighborhoodFactory& nfactory = NeighborhoodFactory::getInstance(ndim);
	vector<vector<int>> nbh = nfactory.neighbor4;
	vector<vector<int>> nbh8 = nfactory.neighbor8;

	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		vector<int> sub = coreParticle2Index(p, ndim);
		for (int n = 0; n < nbh.size(); ++n)
		{
			if (NeighborCheck(sub.begin(), nbh[n].begin(), ndim, dims))
			{
				int idx = Sub2Ind(sub, nbh[n], ndim, dims);
				if (mp[idx] != NULL)
				{
					p->neighbors.push_back(mp[idx]);
				}
			}
		}
		for (int n = 0; n < nbh8.size(); ++n)
		{
			if (NeighborCheck(sub.begin(), nbh8[n].begin(), ndim, dims))
			{
				int idx = Sub2Ind(sub, nbh8[n], ndim, dims);
				if (mp[idx] != NULL)
				{
					p->neighbors8.push_back(mp[idx]);
				}
			}
		}
	}
	return particles;
}

/*
Cluster particles using disjoint set.
*/
vector<vector<CoreParticle*>>
clusterParticles(vector<CoreParticle*>& particles)
{
	set<CoreParticle*> S;
	vector<TK::Node<CoreParticle*>*> nodes;
	map<CoreParticle*, int> imap;
	for (int i = 0; i < particles.size(); ++i)
	{
		TK::Node<CoreParticle*>* n = TK::makeset(particles[i]);
		nodes.push_back(n);
		imap[particles[i]] = i;
		S.insert(particles[i]);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		for (int j = 0; j < p->neighbors8.size(); ++j)
		{
			CoreParticle* q = p->neighbors8[j];
			if (S.find(q) != S.end())
			{
				TK::merge(nodes[i], nodes[imap[q]]);
			}
		}
	}
	vector<TK::Node<CoreParticle*>*> rep = TK::clusters(nodes);
	vector<vector<CoreParticle*>> group(rep.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = distance(rep.begin(), find(rep.begin(), rep.end(), findset(nodes[i])));
		group[k].push_back(nodes[i]->key);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return group;
}


/*
Based on the L2 distance value, find the neighbor that is most straight.

This version uses 8-neighbors.
*/
CoreParticle*
representativeNeighbor(CoreParticle* p, int gen)
{
	float minerr = std::numeric_limits<float>::infinity();
	CoreParticle* best = NULL;
	for (int i = 0; i < p->neighbors8.size(); ++i)
	{
		CoreParticle* q = p->neighbors8[i];
		float len = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
		if (q->gen == gen)
		{
			float err = Abs(p->dval + len - q->dval);
			if (err < minerr)
			{
				minerr = err;
				best = q;
			}
		}
	}
	return best;
}

vector<CoreParticle*>
selectNeighbors(CoreParticle*p, int gen, bool bEight = true)
{
	vector<CoreParticle*> vp;
	if (bEight)
	{
		for (int i = 0; i < p->neighbors8.size(); ++i)
		{
			if (p->neighbors8[i]->gen == gen)
			{
				vp.push_back(p->neighbors8[i]);
			}
		}
	}
	else
	{
		for (int i = 0; i < p->neighbors.size(); ++i)
		{
			if (p->neighbors[i]->gen == gen)
			{
				vp.push_back(p->neighbors[i]);
			}
		}
	}
	return vp;
}

vector<float>
centroid(vector<CoreParticle*>& vp)
{
	vector<float> v(4, 0.0f);
	for (int j = 0; j < vp.size(); ++j)
	{
		v[0] += vp[j]->x;
		v[1] += vp[j]->y;
		v[2] += vp[j]->z;
		v[3] += vp[j]->t;
	}
	v[0] /= vp.size(); v[1] /= vp.size(); v[2] /= vp.size(); v[3] /= vp.size();
	return v;
}

/*
find a particle among VP that is closest to C with 4 elements.
*/
CoreParticle*
closestParticle(vector<CoreParticle*>& vp, vector<float> c)
{
	//vector<float> v = centroid(vp);
	CoreParticle* rep = NULL;
	float mind = std::numeric_limits<float>::infinity();
	for (int j = 0; j < vp.size(); ++j)
	{
		float d = length(c[0] - vp[j]->x, c[1] - vp[j]->y, c[2] - vp[j]->z, c[3] - vp[j]->t);
		if (d < mind)
		{
			mind = d;
			rep = vp[j];
		}
	}
	return rep;
}

CoreParticle*
descendUntilCore(CoreParticle* p)
{
	vector<CoreParticle*> Q(1, p);
	while (Q.empty() == false)
	{
		set<CoreParticle*> S;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* q = Q[i];
			if (q->descendents.empty()) return q;

			for (set<CoreParticle*>::iterator it = q->descendents.begin(); it != q->descendents.end(); it++)
			{
				S.insert(*it);
			}
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
	return NULL;
}

vector<CoreParticle*>
closestDescendents(CoreParticle* p)
{
	vector<CoreParticle*> vn = selectNeighbors(p, p->gen + 1);
	vector<CoreParticle*> vq;
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < vn.size(); ++i)
	{
		CoreParticle* q = vn[i];
		if (q->gen == p->gen && q->id < p->id) continue;

		float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
		if (d < mind)
		{
			mind = d;
			vq.clear();
			vq.push_back(q);
		}
		else if (Abs(d - mind) < 0.001f)
		{
			vq.push_back(q);
		}
	}
	return vq;
}

vector<CoreParticle*>
closestInSameGeneration(CoreParticle* p)
{
	vector<CoreParticle*> vq;
	vector<CoreParticle*> vn = selectNeighbors(p, p->gen);
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < vn.size(); ++i)
	{
		CoreParticle* q = vn[i];
		//if (q->gen == p->gen && q->id < p->id) continue;
		for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
		{
			CoreParticle* r = *it;
			float d = length(2 * p->x - r->x - q->x, 2 * p->y - r->y - q->y, 2 * p->z - r->z - q->z, 2 * p->t - r->t - q->t);
			if (d < mind)
			{
				mind = d;
				vq.clear();
				vq.push_back(q);
			}
			else if (Abs(d - mind) < 0.001f)
			{
				vq.push_back(q);
			}
		}
	}
	for (int i = vq.size() - 1; i >= 0; i--)
	{
		if (vq[i]->id < p->id)
		{
			vq.erase(vq.begin() + i);
		}
	}
	return vq;
}

/*
Link particles without descendents to its neighbors with the same generation that is strong medial particle.
This is to link particles inside isolated local maximum components.
*/
void
propagateDescendency2(vector<CoreParticle*>& P, int ndim)
{
	vector<CoreParticle*> U; //particles that can links to decendency
	for (int i = 0; i < P.size(); ++i)
	{
		if (P[i]->descendents.empty())
		{
			U.push_back(P[i]);
		}
	}
	vector<vector<CoreParticle*>> C = clusterParticles(U);
	vector<CoreParticle*> R;
	for (int i = 0; i < C.size(); ++i)
	{
		CoreParticle* rep = closestParticle(C[i], centroid(C[i]));
		if (rep)
		{
			rep->descendents.insert(rep); //temporary add itself as a descendent
			R.push_back(rep);
		}
	}
	vector<CoreParticle*> Q = R;
	while (Q.empty() == false)
	{
		set<CoreParticle*> S; //particles that are to be lined to decendency
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			for (int j = 0; j < p->neighbors8.size(); ++j)
			{
				CoreParticle* q = p->neighbors8[j];
				if (p->gen == q->gen)
				{
					if (q->descendents.size() == 0)
					{
						S.insert(q);
					}
				}
			}
		}
		//for each particle in S, find the best one according to representative neighbor
		for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			CoreParticle* p = *it;
			vector<CoreParticle*> vp;
			vector<CoreParticle*> vn = selectNeighbors(p, p->gen);
			float mind = std::numeric_limits<float>::infinity();
			for (int j = 0; j < vn.size(); ++j)
			{
				CoreParticle* q = vn[j];
				if (q->descendents.empty() == false)
				{
					float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
					if (d < mind)
					{
						vp.clear();
						vp.push_back(q);
						mind = d;
					}
					else if (Abs(d - mind) < 0.001f)
					{
						vp.push_back(q);
					}
				}
			}
			for (int j = 0; j < vp.size(); ++j)
			{
				CoreParticle* q = vp[j];
				p->descendents.insert(q);
				q->ascendents.insert(p);
			}
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
	for (int i = 0; i < R.size(); ++i)
	{
		R[i]->descendents.clear();
	}
}

/*
Link particles without descendents to its neighbors with the same generation and with a decendent.
*/
void
propagateDescendency(vector<CoreParticle*>& P)
{
	vector<CoreParticle*> Q; //particles that can links to decendency
	for (int i = 0; i < P.size(); ++i)
	{
		if (P[i]->descendents.size()>0)
		{
			Q.push_back(P[i]);
		}
	}
	while (Q.empty() == false)
	{
		set<CoreParticle*> S; //particles that are to be linked to decendency
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			for (int j = 0; j < p->neighbors8.size(); ++j)
			{
				CoreParticle* q = p->neighbors8[j];
				if (p->gen == q->gen)
				{
					if (q->descendents.size() == 0)
					{
						S.insert(q);
					}
				}
			}
		}
		//for each particle in S, find the best one according to representative neighbor
		for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			CoreParticle* p = *it;
			vector<CoreParticle*> vp;
			vector<CoreParticle*> vn = selectNeighbors(p, p->gen);
			float mind = std::numeric_limits<float>::infinity();
			for (int j = 0; j < vn.size(); ++j)
			{
				CoreParticle* q = vn[j];
				if (q->descendents.empty() == false)
				{
					float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
					if (d < mind)
					{
						vp.clear();
						vp.push_back(q);
						mind = d;
					}
					else if (Abs(d - mind) < 0.001f)
					{
						vp.push_back(q);
					}
				}
			}
			for (int j = 0; j < vp.size(); ++j)
			{
				CoreParticle* q = vp[j];
				p->descendents.insert(q);
				q->ascendents.insert(p);
			}
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
}

/*
from each surface particle, move towrad the center and establish ascendent/descendent relation.
*/
vector<CoreParticle*>
propagateParticles(
	vector<CoreParticle*>& particles,
	vector<CoreParticle*>& mp,
	vector<int>& S,
	int ndim,
	const int* dims)
{
	vector<CoreParticle*> Q;
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		p->value = 0;
		if (surfaceParticle(p, ndim, dims))
		{
			Q.push_back(p);
		}
	}
	printf("There are %d surface voxels.\n", Q.size());
	int gen = 1;
	while (Q.empty() == false)
	{
		//vector<CoreParticle*> core;
		for (int i = 0; i < Q.size(); ++i)
		{
			SetVoxel(S, Q[i], gen, ndim, dims);
			Q[i]->gen = gen;
			Q[i]->value = 2;
		}

		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			CoreParticle* p2 = NULL;
			for (int n = 0; n < p->neighbors.size(); ++n)
			{
				CoreParticle* q = p->neighbors[n];
				int sval = GetVoxel(S, q, 0, ndim, dims);
				if (sval == 0)
				{
					Q2.insert(q);
					q->value = 1;
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());

		gen++;
	}
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		vector<CoreParticle*> vq = closestDescendents(p);
		for (int j = 0; j < vq.size(); ++j)
		{
			p->descendents.insert(vq[j]);
			vq[j]->ascendents.insert(p);
		}
	}
	/*map<CoreParticle*, vector<CoreParticle*>> msp;
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		if (p->descendents.empty())
		{
			vector<CoreParticle*> vq = closestInSameGeneration(p);
			if (vq.empty() == false)
			{
				msp[p] = vq;
			}
		}
	}
	for (map<CoreParticle*, vector<CoreParticle*>>::iterator it = msp.begin(); it != msp.end(); ++it)
	{
		CoreParticle* p = (it)->first;
		vector<CoreParticle*> vq = it->second;
		for (int j = 0; j < vq.size(); ++j)
		{
			p->descendents.insert(vq[j]);
			vq[j]->ascendents.insert(p);
		}
	}*/

	//propagateDescendency(particles);
	//propagateDescendency2(particles, ndim);

	//remove non-essential link
	/*while (prune)
	{
	bool bChanged = false;
	for (int i = 0; i < particles.size(); ++i)
	{
	CoreParticle* p = particles[i];
	if (p->ascendents.size()>1)
	{
	for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
	{
	CoreParticle* q = *it;
	if (q->descendents.size() > 1)
	{
	p->ascendents.erase(find(p->ascendents.begin(), p->ascendents.end(), q));
	q->descendents.erase(find(q->descendents.begin(), q->descendents.end(), p));
	bChanged = true;
	break;
	}
	}
	}
	}
	if (!bChanged) break;
	}*/

	return particles;
}

vector<TK::Vertex<CoreParticle*>*>
makeGraphStructure(vector<CoreParticle*>& mp, vector<int>& S,
	int ndim, const int* dims)
{
	vector<CoreParticle*> core;
	vector<CoreParticle*> particles;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i])
		{
			mp[i]->label = 0;
			mp[i]->value = 0;
			particles.push_back(mp[i]);
		}
	}
	labelStrongCoreParticles(particles, ndim, dims);

	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i]->core > 0)
		{
			core.push_back(particles[i]);
			particles[i]->src = particles[i];
		}
		else {
			particles[i]->src = NULL;
		}
	}

	propagateDescendency(particles);
	propagateDescendency2(particles, ndim);

	set<int> sgen;
	for (int i = 0; i < core.size(); ++i)
	{
		sgen.insert(core[i]->gen);
		core[i]->label = i + 1;
	}
	vector<int> vgen;
	vgen.insert(vgen.begin(), sgen.begin(), sgen.end());
	sort(vgen.begin(), vgen.end());
	for (int ig = 0; ig < vgen.size(); ++ig)
	{
		int gval = vgen[ig];
		vector<CoreParticle*> Q;
		for (int i = 0; i < core.size(); ++i)
		{
			if (core[i]->gen == gval)
			{
				Q.push_back(core[i]);
			}
		}
		while (Q.empty() == false)
		{
			set<CoreParticle*> Q2;
			for (int i = 0; i < Q.size(); ++i)
			{
				CoreParticle* p = Q[i];
				SetVoxel(S, p, p->label, ndim, dims);
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
					//for (vector<CoreParticle*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->label <= 0)
					{
						q->label = p->label;
						q->value = p->value + 1;
						q->src = p->src;
						Q2.insert(q);
					}
				}
			}
			Q.clear();
			Q.insert(Q.end(), Q2.begin(), Q2.end());
		}
	}
	//now find labels that are adjacent to each other
	TK::GraphFactory<CoreParticle*>& factory = TK::GraphFactory<CoreParticle*>::GetInstance();
	vector<TK::Vertex<CoreParticle*>*> vertices;
	for (int i = 0; i < core.size(); ++i)
	{
		TK::Vertex<CoreParticle*>* u = factory.makeVertex(core[i]);
		vertices.push_back(u);
		core[i]->vertex = u;
	}
	vector<TK::Edge<CoreParticle*>*> edges;
	//map<pair<CoreParticle*, CoreParticle*>, Edge<CoreParticle*>*> emap;
	int edge_count = 0;
	for (int i = 0; i < mp.size(); ++i)
	{
		CoreParticle* p = mp[i];
		if (p)
		{
			for (int j = 0; j < p->neighbors.size(); ++j)
			{
				CoreParticle* q = p->neighbors[j];
				if ((q->x + 1 == 53 && q->y + 1 == 17 && q->z + 1 == 30) || (p->x + 1 == 53 && p->y + 1 == 17 && p->z + 1 == 30))
				{
					p->id += 0;
				}
				if (p->label > 0 && q->label > 0 && p->label != q->label)
				{
					TK::Vertex<CoreParticle*>* u = vertices[p->label - 1];
					TK::Vertex<CoreParticle*>* v = vertices[q->label - 1];
					TK::Edge<CoreParticle*>* uv = u->findEdge(v);
					TK::Edge<CoreParticle*>* vu = v->findEdge(u);
					float w = weight(p, q);
					if (uv == NULL)
					{
						uv = factory.makeEdge(u, v, w);
						vu = factory.makeEdge(v, u, w);
						u->Add(uv);
						v->Add(vu);
						edge_count += 2;
						edges.push_back(uv);
						edges.push_back(vu);
					}
					if (w < uv->w)
					{
						uv->w = w;
						vu->w = w;
					}
				}
			}
		}
	}
	//vector<Edge<CoreParticle*>*> mst = Kruskal(vertices);
	//printf("#vertices = %d, #edges = %d\n", vertices.size(), edge_count);

	/*for (int i = 0; i < mst.size(); ++i)
	{
	printf("%d %3.3f -- %d %3.3f ==> %f\n",
	mst[i]->u->key->id, mst[i]->u->key->dval, mst[i]->v->key->id, mst[i]->v->key->dval, mst[i]->w);
	}*/
	return vertices;
}

/*
Calculate the saliency of each particle once and for all.
The saliency is calculated as:
	1.0 for surface particles
	sum of (saliency of ascendent / # of descendents) for non-surface particles
*/
void
evaluateSaliency(vector<CoreParticle*>& particles, int ndim, const int* dims)
{
	set<CoreParticle*> Q;
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		if (surfaceParticle(p, ndim, dims))
		{
			p->saliency = 1.0;
			for (set<CoreParticle*>::iterator it = p->descendents.begin(); it != p->descendents.end(); it++)
			{
				Q.insert(*it);
			}
		}
	}
	while (Q.empty() == false)
	{
		set<CoreParticle*> S;
		for (set<CoreParticle*>::iterator jt=Q.begin(); jt != Q.end(); ++jt)
		{
			CoreParticle* p = *jt;
			float sum = 0;
			for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
			{
				CoreParticle* q = *it;
				float val = q->saliency / q->descendents.size();
				sum += val;
			}
			p->saliency = sum;
			for (set<CoreParticle*>::iterator it = p->descendents.begin(); it != p->descendents.end(); it++)
			{
				S.insert(*it);
			}
		}
		Q.clear();
		for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)

		{
			CoreParticle* q = *it;
			if (q->saliency != q->saliency)
			{
				Q.insert(q);
			}
		}
	}
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		if (p->saliency != p->saliency) //still NaN. Is something wrong?
		{
			p->saliency = 0.0f;
			printf("NaN Saliency: %d(%d,%d), %d, %f, %d, %d\n", 
				p->id, p->x, p->y, p->gen, p->dval, p->ascendents.size(), p->descendents.size());
		}
	}
}

void
partition(vector<TK::Edge<CoreParticle*>*>& tree, vector<CoreParticle*>& particles, vector<int>& S,
	set<TK::Edge<CoreParticle*>*>& toremove, int ndim, const int* dims)
{
	set<TK::Vertex<CoreParticle*>*> vertices;
	for (int i = 0; i < tree.size(); ++i) {
		vertices.insert(tree[i]->u);
		vertices.insert(tree[i]->v);
	}
	vector<TK::Node<TK::Vertex<CoreParticle*>*>*> nodes;
	map<TK::Vertex<CoreParticle*>*, TK::Node<TK::Vertex<CoreParticle*>*>*> vnmap;
	for (set<TK::Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); it++) {
		TK::Node<TK::Vertex<CoreParticle*>*>* n = TK::makeset(*it);
		nodes.push_back(n);
		vnmap[*it] = n;
	}
	for (int i = 0; i < tree.size(); ++i) {
		if (toremove.find(tree[i]) == toremove.end())
		{
			merge(vnmap[tree[i]->u], vnmap[tree[i]->v]);
		}
	}
	vector<TK::Node<TK::Vertex<CoreParticle*>*>*> reps = clusters(nodes);
	map<TK::Node<TK::Vertex<CoreParticle*>*>*, int> nimap;
	for (int i = 0; i < reps.size(); ++i)
	{
		nimap[reps[i]] = i + 1;
	}
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
	}
	set<int> sgen;
	for (set<TK::Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); ++it)
	{
		TK::Node<TK::Vertex<CoreParticle*>*>* n = vnmap[*it];
		TK::Node<TK::Vertex<CoreParticle*>*>* r = findset(n);
		n->key->key->label = nimap[r];
		sgen.insert(n->key->key->gen);
	}
	vector<int> vgen;
	vgen.insert(vgen.begin(), sgen.begin(), sgen.end());
	sort(vgen.begin(), vgen.end());
	for (int ig = 0; ig < vgen.size(); ++ig)
	{
		int gval = vgen[ig];
		vector<CoreParticle*> Q;
		for (set<TK::Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); ++it)
		{
			if ((*it)->key->gen == gval)
			{
				Q.push_back((*it)->key);
			}
		}
		while (Q.empty() == false)
		{
			set<CoreParticle*> Q2;
			for (int i = 0; i < Q.size(); ++i)
			{
				CoreParticle* p = Q[i];
				SetVoxel(S, p, p->label, ndim, dims);
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->label <= 0)
					{
						q->label = p->label;
						q->value = p->value + 1;
						Q2.insert(q);
					}
				}
			}
			Q.clear();
			Q.insert(Q.end(), Q2.begin(), Q2.end());
		}
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "GroupPoints", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [A B C D] = MovingDotGrouping(P)");
		return;
	}

	vector<unsigned char> L;
	int ndimL;
	const int* dimsL;
	{
		mxClassID classIdP;
		LoadData(L, prhs[0], classIdP, ndimL, &dimsL);
	}
	float thres = 900.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[1], classMode);
	}
	printf("Threshold = %f\n", thres);
	int nvoxels = numberOfElements(ndimL, dimsL);
	vector<CoreParticle*> mp = generateParticleMap(L, ndimL, dimsL);
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndimL, dimsL);
	vector<float> D(L.size(), 0.0f);
	{
		vector<unsigned char> iL(L.size(), 0);
		for (int i = 0; i < L.size(); ++i)
		{
			iL[i] = L[i] ? 0 : 1;
		}
		vector<float> vs(ndimL, 1.0f);
		DistanceTransformEuclidF(D, iL, vs, ndimL, dimsL);
		for (int i = 0; i < D.size(); ++i)
		{
			if (mp[i] != NULL)
			{
				mp[i]->dval = D[i];
			}
		}
		iL.clear(); //to save memory 
		D.clear(); //to save memory 
	}
	vector<int> S(nvoxels, 0);
	propagateParticles(particles, mp, S, ndimL, dimsL);
	evaluateSaliency(particles, ndimL, dimsL);
	vector<int> S2(nvoxels, 0);
	vector<TK::Vertex<CoreParticle*>*> vertices = makeGraphStructure(mp, S2, ndimL, dimsL);
	vector<TK::Edge<CoreParticle*>*> mst = Kruskal(vertices);

	set<TK::Edge<CoreParticle*>*> toremove;
	for (int i = 0; i < mst.size(); ++i)
	{
		if (mst[i]->w > thres)
		{
			toremove.insert(mst[i]);
		}
	}
	vector<int> S3(nvoxels, 0);
	partition(mst, particles, S3, toremove, ndimL, dimsL);

	if (nlhs >= 1)
	{
		int dims[] = { particles.size(), 12 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < particles.size(); ++i)
		{
			CoreParticle* p = particles[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->z);
			SetData2(F, i, 3, dims[0], dims[1], p->t);
			SetData2(F, i, 4, dims[0], dims[1], p->label);
			SetData2(F, i, 5, dims[0], dims[1], p->gen);
			SetData2(F, i, 6, dims[0], dims[1], p->id);
			SetData2(F, i, 7, dims[0], dims[1], surfaceParticle(p, ndimL, dimsL) ? 1 : 0);
			SetData2(F, i, 8, dims[0], dims[1], p->core);
			SetData2(F, i, 9, dims[0], dims[1], (int)(p->saliency * 100));
			SetData2(F, i, 10, dims[0], dims[1], (int)p->descendents.size());
			SetData2(F, i, 11, dims[0], dims[1], (int)p->ascendents.size());
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		vector<TK::Edge<CoreParticle*>*> all = GetAllEdges(vertices);
		vector<TK::Edge<CoreParticle*>*> edges;
		for (int i = 0; i < all.size(); ++i)
		{
			if (all[i]->u->key->id < all[i]->v->key->id)
			{
				edges.push_back(all[i]);
			}
		}
		int dims[] = { edges.size(), 4 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < edges.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], edges[i]->u->key->id);
			SetData2(F, i, 1, dims[0], dims[1], edges[i]->v->key->id);
			SetData2(F, i, 2, dims[0], dims[1], (int)(100 * edges[i]->w));
			SetData2(F, i, 3, dims[0], dims[1], (int)(edges[i]->type));
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreData(S, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreData(S2, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 5)
	{
		int dims[] = { mst.size(), 4 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < mst.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], mst[i]->u->key->id);
			SetData2(F, i, 1, dims[0], dims[1], mst[i]->v->key->id);
			SetData2(F, i, 2, dims[0], dims[1], (int)(mst[i]->w));
			SetData2(F, i, 3, dims[0], dims[1], (int)(mst[i]->type));
		}
		plhs[4] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 6)
	{
		plhs[5] = StoreData(S3, mxINT32_CLASS, ndimL, dimsL);
	}
	/*if (nlhs >= 5)
	{
		vector<pair<int, int>> pairs;
		for (int i = 0; i < particles.size(); ++i)
		{
			CoreParticle* p = particles[i];
			//if (surfaceParticle(p, ndimL) == false) continue; //keep it only for surface particles
			for (set<CoreParticle*>::iterator it = p->descendents.begin(); it != p->descendents.end(); ++it)
			{
				pairs.push_back(pair<int, int>(p->id, (*it)->id));
			}
		}
		int dims[] = { pairs.size(), 2 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < pairs.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], pairs[i].first);
			SetData2(F, i, 1, dims[0], dims[1], pairs[i].second);
		}
		plhs[4] = StoreData(F, mxINT32_CLASS, 2, dims);
	}*/
	TK::GraphFactory<CoreParticle*>::GetInstance().Clean();
	CoreParticleFactory::getInstance().clean();

	mexUnlock();
}

