#include <CoreParticleMakeGraph.h>
#include <NeighborhoodFactory.h>
#include <GraphFactory.h>
#include <szMiscOperations.h>

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

vector<CoreParticle*>
closestDescendents(CoreParticle* p, bool bEightNeighbor)
{
	vector<CoreParticle*> vn = selectNeighbors(p, p->gen + 1, bEightNeighbor);
	vector<CoreParticle*> vq;
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < vn.size(); ++i)
	{
		CoreParticle* q = vn[i];
		//if (q->gen == p->gen && q->id < p->id) continue;

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
closestAscendents(CoreParticle* p, bool bEightNeighbor)
{
	vector<CoreParticle*> vn = selectNeighbors(p, p->gen - 1, bEightNeighbor);
	vector<CoreParticle*> vq;
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < vn.size(); ++i)
	{
		CoreParticle* q = vn[i];
		//if (q->gen == p->gen && q->id < p->id) continue;

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

	return particles;
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
			mp[i]->vertex = NULL;
			particles.push_back(mp[i]);
		}
	}
	//labelStrongCoreParticles(particles, ndim, dims);

	for (int i = 0; i < particles.size(); ++i)
	{
		if (strongMedialParticle(particles[i], ndim, dims) > 0)
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
	return vertices;
}

