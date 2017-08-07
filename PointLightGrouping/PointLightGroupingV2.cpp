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
#include <MCpoint.h>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>
#include <PointLight.h>
#include <DisjointSet.h>
#include <Hungarian.h>

using namespace _PointLightNS;

int PointLight::_id = 0;

vector<int>
clusterPoints(vector<PointLight*>& points)
{
	vector<Node<PointLight*>*> nodes;
	map<PointLight*, int> pmap;
	for (int i = 0; i < points.size(); ++i)
	{
		nodes.push_back(makeset(points[i]));
		pmap[points[i]] = i;
	}
	for (int i = 0; i < points.size(); ++i)
	{
		PointLight::PointLightTriple w = points[i]->winner();
		if (w.p != NULL)
		{
			PointLight::PointLightTriple wp = w.p->winner();
			if (wp.r == points[i]) 
			{
				merge(nodes[i], nodes[pmap[w.p]]);
			}
		}
		if (w.r != NULL)
		{
			PointLight::PointLightTriple wr = w.r->winner();
			if (wr.p == points[i])
			{
				merge(nodes[i], nodes[pmap[w.r]]);
			}
		}
	}
	vector<Node<PointLight*>*> reps = clusters(nodes);
	map<PointLight*, int> cmap;
	for (int i = 0; i < reps.size(); ++i)
	{
		cmap[reps[i]->key] = i;
	}
	vector<int> labels(points.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = cmap[findset(nodes[i])->key];
		labels[i] = k;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		PointLight* p = nodes[i]->key;
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			PointLight* q = nodes[j]->key;
		}
		int k = cmap[findset(nodes[i])->key];
		labels[i] = k;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return labels;
}

vector<vector<PointLight*>>
organizeIntoFrames(vector<PointLight*>& P)
{
	vector<pair<float, PointLight*>> pairs(P.size());
	for (int i = 0; i < P.size(); ++i)
	{
		pairs[i].first = P[i]->frame;
		pairs[i].second = P[i];
	}
	sort(pairs.begin(), pairs.end());
	vector<vector<PointLight*>> frames;
	vector<PointLight*> frame;
	for (int i = 0; i < pairs.size(); ++i)
	{
		if (frame.empty() || frame[0]->frame == pairs[i].first)
		{
			frame.push_back(pairs[i].second);
		}
		else
		{
			frames.push_back(frame);
			frame.clear();
			frame.push_back(pairs[i].second);
		}
	}
	frames.push_back(frame);
	return frames;
}

/*
Check if there is already a point that is closely aligned to the line p-q and is temporarily closer to p.
*/
bool
isRedundant(PointLight* p, PointLight* q, vector<PointLight*>& Q, float thres)
{
	CParticleF p0(p->x, p->y);
	CParticleF q0(q->x, q->y);
	for (int i = 0; i < Q.size(); ++i)
	{
		CParticleF x(Q[i]->x, Q[i]->y);
		CParticleF c = Closest2Line(p0, q0, x);
		if (Distance(x, c) < thres)
		{
			return true;
		}
	}
	return false;
}

float
linkMeasure(PointLight* q)
{
	float measure = 0;
	while (q != NULL)
	{
		PointLight* next = NULL;
		PointLight::PointLightTriple tp = q->winner();
		if (tp.p != NULL && tp.p != q)
		{
			PointLight::PointLightTriple tp2 = tp.p->winner();
			if (tp2.r == q)
			{
				measure += 1.0f;
				next = tp.p;
			}
		}
		q = next;
	}
	return measure;
}

void
addFrame(vector<PointLight*>& Q, vector<PointLight*>& frame, int range, float thres, float sigma)
{
	vector<PointLight*> P;
	vector<vector<float>> C0;
	vector<vector<int>> cindx;
	for (int j = 0; j < Q.size(); ++j)
	{
		PointLight* p = Q[j];
		if (p->state != Connected && p->state != Appeared) //need an extension
		{
			vector<float> vals(frame.size(), 0);
			vector<int> ids(frame.size(), -1);
			bool bKeep = false;
			for (int i = 0; i < frame.size(); ++i)
			{
				PointLight* q = frame[i];
				if (p->frame < q->frame - range) continue;

				if (p->id == 559 && q->id == 593)
				{
					i += 0;
				}
				for (int k = 0; k < p->candidates.size(); ++k)
				{
					PointLight::Candidate c(p, q, q);
					float val = 0;
					float fit = 0;
					if (PointLight::isStationary(p->candidates[k]))
					{
						fit = PointLight::_CompatibilityStationary(p->candidates[k], c);
						val = fit * p->candidates[k].prob;
					}
					else
					{
						PointLight::Candidate c2(p->candidates[k].tp.p, p->candidates[k].tp.q, q);
						float lnk = linkMeasure(p);
						fit = PointLight::_Compatibility(c2, c);
						val = (fit + lnk) * p->candidates[k].prob;
					}
					if (val > vals[i])
					{
						ids[i] = k;
						vals[i] = val;
					}
					if (fit > _PointLightNS::Threshold)
					{
						bKeep = true;
					}
				}
			}
			if (bKeep)
			{
				P.push_back(p);
				C0.push_back(vals);
				cindx.push_back(ids);
			}
		}
	}
	printf("dim: %d, %d\n", frame.size(), P.size());
	int n = Max(frame.size(), P.size());
	const int dimsC[] = { n, n };
	//convert cost matrix to a flat vector with appropriate 0-padding
	vector<float> C(n*n, 0.0f);
	for (int i = 0; i < frame.size(); ++i)
	{
		for (int j = 0; j < P.size(); ++j)
		{
			C[i*n + j] = C0[j][i];
		}
	}
	Hungarian algo(C, dimsC);
	float cost = algo.solve();
	

	/*for (int i = 0; i < frame.size(); ++i)
	{
		int j = algo.getJob(i);
		printf("%d=>%d ", frame[i]->id, j < P.size() ? P[j]->id : -1);
	}
	printf("\n");*/

	for (int i = 0; i < frame.size(); ++i)
	{
		PointLight* q = frame[i];
		int j = algo.getJob(i);
		if (j < P.size()) // && algo.getCost(j, i) > _PointLightNS::Threshold)
		{
			PointLight* p = P[j];

			q->candidates0.push_back(PointLight::Candidate(p, q, q));
			int k = cindx[j][i];
			p->candidates0.push_back(PointLight::Candidate(p->candidates[k].tp.p, p->candidates[k].tp.q, q));
			/*printf(">> (%d,%d,%d), (%d,%d,%d)\n",
				p->id, q->id, q->id, p->candidates[k].tp.p->id, p->candidates[k].tp.q->id, q->id);*/
		}
		else
		{
			q->candidates0.push_back(PointLight::Candidate(q, q, q));
			//printf(">> (%d,%d,%d)\n", q->id, q->id, q->id);
		}
	}
	Q.insert(Q.end(), frame.begin(), frame.end());
}

void
update(vector<PointLight*>& P, int niter)
{
	for (int i = 0; i < P.size(); ++i)
	{
		P[i]->updateCandidates();
	}
	for (int n = 0; n < niter; ++n)
	{
		for (int i = 0; i < P.size(); ++i)
		{
			P[i]->updateFitness();
			P[i]->updateProb();
		}
		for (int i = 0; i < P.size(); ++i)
		{
			//P[i]->updateProb();
		}
	}
	for (int i = 0; i < P.size(); ++i)
	{
		P[i]->updateState();
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "GroupPoints", __DATE__, __TIME__);
	if (nrhs < 0 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [A B C D] = GroupPoints(numiter, noise_var)");
		return;
	}
	PointLight::_id = 0;

	vector<PointLight*> P;
	{
		vector<float> P0;
		const int* dimsP;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			float fr = GetData2(P0, i, 2, dimsP[0], dimsP[1], (float)0);
			float gid = GetData2(P0, i, 3, dimsP[0], dimsP[1], (float)0);
			PointLight* pl = new PointLight(x, y, fr, gid);
			P.push_back(pl);
		}
	}

	int numIter = 1;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[1], classMode);
	}
	float sigma = 10.0;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(sigma, prhs[2], classMode);
	}
	float rate = 0.25;
	vector<vector<PointLight*>> frames = organizeIntoFrames(P);
	int numFrames = frames.size();
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(numFrames, prhs[3], classMode);
	}

	vector<PointLight*> Q;
	for (int k = 0; k < numFrames; ++k)
	{
		printf("Frame: %d\n", k);
		addFrame(Q, frames[k], 5, 0.1, sigma);
		update(Q, numIter);
	}
	vector<int> labels = clusterPoints(Q);

	if (nlhs >= 1)
	{
		const int dims[] = { Q.size(), 5 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			PointLight* p = Q[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->frame);
			SetData2(F, i, 3, dims[0], dims[1], (float)p->id);
			SetData2(F, i, 4, dims[0], dims[1], (float)p->gid);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { labels.size(), 1 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], labels[i]);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		int nc = MaxNumCandidates;
		const int dims[] = { Q.size(), nc};
		vector<int> F(dims[0] * dims[1], -1);
		for (int i = 0; i < dims[0]; ++i)
		{
			PointLight* p = Q[i];
			for (int j = 0, k=0; j < P[i]->candidates.size(); ++j)
			{
				SetData2(F, i, k++, dims[0], dims[1], P[i]->candidates[j].tp.p->id);
				SetData2(F, i, k++, dims[0], dims[1], P[i]->candidates[j].tp.q->id);
				SetData2(F, i, k++, dims[0], dims[1], P[i]->candidates[j].tp.r->id);
			}
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	for (int i = 0; i < P.size(); ++i)
	{
		delete P[i];
	}
	mexUnlock();
}

