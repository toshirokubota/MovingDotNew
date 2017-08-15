#pragma once
#include <CoreParticle.h>
#include <CoreParticleUtility.h>
#include <Graph.h>

vector<float>
centroid(vector<CoreParticle*>& vp);

vector<CoreParticle*>
propagateParticles(
	vector<CoreParticle*>& particles,
	vector<CoreParticle*>& mp,
	vector<int>& S,
	int ndim,
	const int* dims);

vector<TK::Vertex<CoreParticle*>*>
makeGraphStructure(vector<CoreParticle*>& mp, vector<int>& S,
	int ndim, const int* dims);
