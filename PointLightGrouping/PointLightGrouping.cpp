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

void
addFrame(vector<PointLight*>& Q, vector<PointLight*>& frame, int range, float thres, float sigma)
{
	for (int i = 0; i < frame.size(); ++i)
	{
		PointLight* q = frame[i];
		vector<pair<float, pair<PointLight*, int>>> pairs;
		vector<pair<float, pair<PointLight*, int>>> pairs2;
		for (int j = 0; j < Q.size(); ++j)
		{
			if (Abs(q->frame - Q[j]->frame) <= range)
			{
				PointLight* p = Q[j];
				if (p->state != Connected && p->state != Appeared) //need an extension
				{
					for (int k = 0; k < p->candidates.size(); ++k)
					{
						PointLight::Candidate c(p, q, q);
						if (PointLight::isStationary(p->candidates[k]))
						{
							float cmp = PointLight::_CompatibilityStationary(p->candidates[k], c) * p->candidates[k].prob;
							pairs2.push_back(pair<float, pair<PointLight*, int>>(cmp, pair<PointLight*, int>(p, k)));
							if ((q->id == 415 || q->id == 416) && (p->id == 404 || p->id == 403))
							{
								printf("addFrameS: %d %d:  %f, %f, %f\n", q->id, p->id,
									PointLight::_CompatibilityStationary(p->candidates[k], c), p->candidates[k].prob, cmp);
							}
						}
						else
						{
							PointLight::Candidate c2(p->candidates[k].tp.p, p->candidates[k].tp.q, q);
							float cmp = PointLight::_Compatibility(c2, c)* p->candidates[k].prob;
							pairs.push_back(pair<float, pair<PointLight*, int>>(cmp, pair<PointLight*, int>(p, k)));
							if ((q->id == 415 || q->id == 416) && (p->id == 404 || p->id==403))
							{
								printf("addFrame: %d %d:  %f, %f, %f\n", q->id, p->id,
									PointLight::_Compatibility(c2, c), p->candidates[k].prob, cmp);
							}
						}
					}
				}
			}
		}
		sort(pairs.begin(), pairs.end());
		sort(pairs2.begin(), pairs2.end());
		for (int j = pairs.size() - 1; j >= Max(0, (int)pairs.size() - MaxNumNeighbors); j--)
		{
			if (pairs[j].first <= thres) break;
			PointLight* p = pairs[j].second.first;
			int k = pairs[j].second.second;
			q->candidates0.push_back(PointLight::Candidate(p, q, q));
			//p->candidates[k].tp.r = q;
			p->candidates0.push_back(PointLight::Candidate(p->candidates[k].tp.p, p->candidates[k].tp.q, q));
		}
		float sthres = 0;
		for (int j = pairs2.size() - 1; j >= 0 && q->candidates.size() < MaxNumNeighbors; j--)
		{
			if (pairs2[j].first <= sthres) break;
			PointLight* p = pairs2[j].second.first;
			int k = pairs2[j].second.second;
			q->candidates0.push_back(PointLight::Candidate(p, q, q));

			p->candidates0.push_back(PointLight::Candidate(p->candidates[k].tp.p, p->candidates[k].tp.q, q));
			if (j == pairs2.size() - 1)
			{
				sthres = pairs2[j].first / 2.0f;
			}
		}		
		q->candidates0.push_back(PointLight::Candidate(q, q, q)); //always add an isolated candidate
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
	vector<int> labels = clusterPoints(P);

	if (nlhs >= 1)
	{
		const int dims[] = { P.size(), 5 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			PointLight* p = P[i];
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
		const int dims[] = { P.size(), nc};
		vector<int> F(dims[0] * dims[1], -1);
		for (int i = 0; i < dims[0]; ++i)
		{
			PointLight* p = P[i];
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

