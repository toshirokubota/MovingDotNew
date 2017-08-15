#include <CoreParticleUtility.h>
#include <szMexUtility.h>
#include <NeighborhoodFactory.h>
#include <DisjointSet.h>
#include <szDistanceTransformNonIsotropic.h>
#include <map>
using namespace std;

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
		//mexErrMsgTxt("particleNdim: unsupported number of dimensions. It has to be between 1 and 4.");
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
		//mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
	}
	return idx;
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
	{
		vector<float> D(L.size(), 0.0f);
		vector<unsigned char> iL(L.size(), 0);
		vector<float> vs(ndim, 1.0f);
		DistanceTransformEuclidF(D, iL, vs, ndim, dims);
		for (int i = 0; i < D.size(); ++i)
		{
			if (mp[i] != NULL)
			{
				mp[i]->dval = D[i];
			}
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
				merge(nodes[i], nodes[imap[q]]);
			}
		}
	}
	vector<TK::Node<CoreParticle*>*> rep = clusters(nodes);
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
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1 || p->z == 0 || p->z == dims[2] - 1)
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
	/*if (_EightNeighbor)
	{
		ns = NeighborhoodFactory::getInstance(ndim).neighbor8.size();
	}*/
	return p->neighbors.size() < ns;
}

bool
medialParticle(CoreParticle* p, int ndim)
{
	return p->ascendents.size() >= 2 * (ndim - 1) || p->descendents.size() == 0;
}

/*
Ascendents of p are not linearly separable.
*/
bool
strongMedialParticle(CoreParticle* p, int ndim, const int* dims)
{
	if (p->core >= 0) return p->core>0;
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

bool
inflectionParticle(CoreParticle* p, int ndim)
{
	//return surfaceParticle(p, ndim) && p->descendents.size() >= 2;  //ndim;
	return false; //disabled -- 08/08/2017
}//
