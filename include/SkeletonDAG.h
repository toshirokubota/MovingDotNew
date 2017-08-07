#pragma once
#include <Graph.h>
#include <StationaryParticle.h>
#include <mex.h>
#include <map>

class SkeletonDAG
{ 
public:
	static SkeletonDAG fromMxArray(const mxArray* points, const mxArray* edges);
	mxArray* PointsToMxArray();
	mxArray* EdgesToMxArray();
	bool add(StationaryParticle* p);
	bool connect(StationaryParticle* p, StationaryParticle* q, float  w=0.0f);
	vector<Vertex<StationaryParticle*>*> vertices;
	vector<Edge<StationaryParticle*>*> edges;
	map<StationaryParticle*, int> vmap;
};
