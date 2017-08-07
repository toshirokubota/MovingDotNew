#include <SkeletonDAG.h>
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <GraphFactory.h>

bool
SkeletonDAG::add(StationaryParticle* p)
{
	if (vmap.find(p) != vmap.end())
	{
		return false;
	}
	else
	{
		GraphFactory<StationaryParticle*>& factory = GraphFactory<StationaryParticle*>::GetInstance();
		vertices.push_back(factory.makeVertex(p));
		return true;
	}
}

bool
SkeletonDAG::connect(StationaryParticle* p, StationaryParticle* q, float w)
{
	if (vmap.find(p) == vmap.end() || vmap.find(q) == vmap.end())
	{
		return false;
	}
	else
	{
		GraphFactory<StationaryParticle*>& factory = GraphFactory<StationaryParticle*>::GetInstance();
		int u = vmap[p];
		int v = vmap[q];
		Edge<StationaryParticle*>* edge = factory.makeEdge(vertices[u], vertices[v], w, Forward);
		edges.push_back(edge);
		return true;
	}
}

SkeletonDAG
SkeletonDAG::fromMxArray(const mxArray* pptr, const mxArray* eptr)
{
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	GraphFactory<StationaryParticle*>& gfactory = GraphFactory<StationaryParticle*>::GetInstance();
	SkeletonDAG skeleton;

	const int* dimsP;
	vector<float> P0;
	mxClassID classIdP;
	int ndimP;
	LoadData(P0, pptr, classIdP, ndimP, &dimsP);
	for (int i = 0; i < dimsP[0]; ++i)
	{
		float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
		float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
		StationaryParticle* p = sfactory.makeParticle(CParticleF(x, y));
		skeleton.vertices.push_back(gfactory.makeVertex(p));
		skeleton.vmap[p] = i;
	}
	vector<int> T0;
	mxClassID classIdT;
	int ndimT;
	const int* dimsT;
	LoadData(T0, eptr, classIdT, ndimT, &dimsT);
	for (int i = 0; i < dimsT[0]; ++i)
	{
		int v1 = GetData2(T0, i, 0, dimsT[0], dimsT[1], (int)0) - 1;
		int v2 = GetData2(T0, i, 1, dimsT[0], dimsT[1], (int)0) - 1;
		Edge<StationaryParticle*>* e1 = gfactory.makeEdge(skeleton.vertices[v1], skeleton.vertices[v2]);
		Edge<StationaryParticle*>* e2 = gfactory.makeEdge(skeleton.vertices[v2], skeleton.vertices[v1]);
		skeleton.vertices[v1]->Add(e1);
		skeleton.vertices[v2]->Add(e2);
		skeleton.edges.push_back(e1);
		skeleton.edges.push_back(e2);
	}

	return skeleton;
}


mxArray* SkeletonDAG::PointsToMxArray()
{
	const int dims[] = { vertices.size(), 3 };
	vector<float> F(dims[0] * dims[1]);
	for (int i = 0; i<dims[0]; ++i)
	{
		StationaryParticle* p = vertices[i]->key;
		SetData2(F, i, 0, dims[0], dims[1], p->getP().m_X);
		SetData2(F, i, 1, dims[0], dims[1], p->getP().m_Y);
		SetData2(F, i, 2, dims[0], dims[1], (float)p->getId());
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);

}
mxArray* SkeletonDAG::EdgesToMxArray()
{
	vector<int> idx;
	for (int i = 0; i<edges.size(); ++i)
	{
		StationaryParticle* p = edges[i]->u->key;
		StationaryParticle* q = edges[i]->v->key;
		int pi = vmap[p];
		int qi = vmap[q];
		if (edges[i]->type == Forward || pi < qi)
		{
			idx.push_back(i);
		}
	}
	const int dims[] = { idx.size(), 3 };
	vector<float> F(dims[0] * dims[1]);
	for (int i = 0; i < dims[0]; ++i)
	{
		int j = idx[i];
		StationaryParticle* p = edges[j]->u->key;
		StationaryParticle* q = edges[j]->v->key;
		int pi = vmap[p];
		int qi = vmap[q];
		SetData2(F, i, 0, dims[0], dims[1], (float)(pi + 1));
		SetData2(F, i, 1, dims[0], dims[1], (float)(qi + 1));
		SetData2(F, i, 2, dims[0], dims[1], edges[j]->w);
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}
