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
#include <DisjointSet.h>
#include <Hungarian.h>
#include <SimplePoint.h>

int SimplePoint::_id;

struct SimpleLinkedPoint
{
	SimpleLinkedPoint(SimplePoint* p0 = NULL, SimpleLinkedPoint* next = NULL, SimpleLinkedPoint* prev = NULL)
	{
		p = p0;
		this->next = next;
		this->prev = prev;
		fit = 0;
	}
	SimplePoint* p;
	SimpleLinkedPoint* next;
	SimpleLinkedPoint* prev;
	float fit;
};

vector<int>
clusterPoints(vector<SimpleLinkedPoint*>& points)
{
	vector<Node<SimplePoint*>*> nodes;
	map<SimplePoint*, int> pmap;
	for (int i = 0; i < points.size(); ++i)
	{
		nodes.push_back(makeset(points[i]->p));
		pmap[points[i]->p] = i;
	}
	for (int i = 0; i < points.size(); ++i)
	{
		SimpleLinkedPoint* q = points[i]->next;
		if (q != NULL)
		{
			merge(nodes[i], nodes[pmap[q->p]]);
		}
	}
	vector<Node<SimplePoint*>*> reps = clusters(nodes);
	map<SimplePoint*, int> cmap;
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
		SimplePoint* p = nodes[i]->key;
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			SimplePoint* q = nodes[j]->key;
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

vector<vector<SimpleLinkedPoint*>>
organizeIntoFrames(vector<SimpleLinkedPoint*>& P)
{
	vector<pair<float, SimpleLinkedPoint*>> pairs(P.size());
	for (int i = 0; i < P.size(); ++i)
	{
		pairs[i].first = P[i]->p->frame;
		pairs[i].second = P[i];
	}
	sort(pairs.begin(), pairs.end());
	vector<vector<SimpleLinkedPoint*>> frames;
	vector<SimpleLinkedPoint*> frame;
	for (int i = 0; i < pairs.size(); ++i)
	{
		if (frame.empty() || frame[0]->p->frame == pairs[i].first)
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

float
linkMeasure(SimpleLinkedPoint* q)
{
	float measure = 0;
	int nl = 1;
	float wgt = 1.0f;
	while (q != NULL)
	{
		//measure += q->fit * wgt;
		measure += 1.0 / nl;
		q = q->prev;
		nl++;
		wgt *= 0.95;
	}
	return measure;
}

float
Fitness(SimplePoint* q, SimpleLinkedPoint * p, float sigma)
{
	float x, y, df;
	if (p->prev == NULL)
	{
		x = p->p->x;
		y = p->p->y;
	}
	else if (p->prev->prev == NULL)
	{
		float vx = (p->p->x - p->prev->p->x);
		float vy = (p->p->y - p->prev->p->y);
		x = p->p->x + vx;
		y = p->p->y + vy;
	}
	else
	{
		//fit quadratic splines (at*t+bt+c) using p at t=0, p->prev at t=-1, and p->prev->prev at t=-2,
		//and predicts q at t=1.
		float bx = (p->p->x - p->prev->prev->p->x) / 2.0f;
		float by = (p->p->y - p->prev->prev->p->y) / 2.0f;
		float ax = (p->p->x - 2 * p->prev->p->x + p->prev->prev->p->x) / 2.0f;
		float ay = (p->p->y - 2 * p->prev->p->y + p->prev->prev->p->y) / 2.0f;

		x = p->p->x + bx + ax;
		y = p->p->y + by + ay;
	}
	float d = Distance(q->x, q->y, x, y);
	d = d * d;
	float val = exp(-d / (18.0*sigma*sigma));
	val *= pow(0.9, (float)(q->frame - p->p->frame - 1));
	return val;
}

float
getAdaptiveThreshold(Hungarian& algo)
{
	vector<float> vals;
	for (int i = 0; i < algo.nWorkers; ++i)
	{
		float c = algo.getCost(i, algo.getJob(i));
		if (c > 0)
		{
			vals.push_back(c);
		}
	}
	sort(vals.begin(), vals.end());
	float q1 = vals[vals.size() / 4];
	float q3 = vals[3 * vals.size() / 4];
	float thres = q1 - 1.5 * (q3 - q1);
	return thres;
}

vector<SimpleLinkedPoint*>
addFrame(vector<SimpleLinkedPoint*>& Z, vector<SimpleLinkedPoint*>& frame, float sigma)
{
	vector<SimpleLinkedPoint*> P;
	vector<vector<float>> C0;
	vector<vector<int>> cindx;
	float thres = 1.0e-5;
	for (int j = 0; j < Z.size(); ++j)
	{
		{
			vector<float> vals(frame.size(), 0);
			bool bKeep = false;
			for (int i = 0; i < frame.size(); ++i)
			{
				if (Z[j]->p->id == 157 && frame[i]->p->id == 165)
				{
					j += 0;
				}
				float fit = Fitness(frame[i]->p, Z[j], sigma);
				float val = linkMeasure(frame[i]);
				vals[i] = fit * val;
				if (fit > thres)
				{
					bKeep = true;
				}

			}
			if (bKeep)
			{
				P.push_back(Z[j]);
				C0.push_back(vals);
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
	float sum1 = 0, sum2 = 0;
	for (int i = 0; i < algo.nJobs; ++i)
	{
		printf("%f ", algo.getCost(i, algo.getWorker(i)));
		sum1 += algo.getCost(i, algo.getWorker(i));
	}
	printf("\n");
	for (int i = 0; i < algo.nJobs; ++i)
	{
		printf("%f ", algo.getCost(i, algo.getJob(i)));
		sum2 += algo.getCost(i, algo.getJob(i));
	}
	printf("\n");
	float athres = 0.1; // getAdaptiveThreshold(algo);
	printf("cost=%f, sum1=%f, sum2=%f, athres=%f\n", cost, sum1, sum2, athres);

	for (int i = 0; i < frame.size(); ++i)
	{
		SimpleLinkedPoint* q = frame[i];
		int j = algo.getJob(i);
		if (algo.getCost(i, j) < athres) continue;

		if (j < P.size()) 
		{
			SimpleLinkedPoint* p = P[j];
			p->next = q;
			q->prev = p;
			//p->fit = Fitness(p, sigma);
			q->fit = Fitness(q->p, p, sigma);
			if (q->fit < thres)
			{
				q->fit = 0;
				p->next = NULL;
				q->prev = NULL;
			}
		}
	}
	vector<SimpleLinkedPoint*> Q = frame;
	for (int i = 0; i < Z.size(); ++i)
	{
		if (Z[i]->next == NULL)
		{
			Q.push_back(Z[i]);
		}
	}
	return Q;
}

void
updateSigma(vector<SimpleLinkedPoint*>& P, float& sigma, int& n)
{
	float sum = sigma * n;
	for (int i = 0; i < P.size(); ++i)
	{
		if (P[i]->prev != NULL)
		{
			sum += Distance(P[i]->p->x, P[i]->p->y, P[i]->prev->p->x, P[i]->prev->p->y);
			n++;
		}
	}
	sigma = sum / n;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "GroupPoints", __DATE__, __TIME__);
	if (nrhs < 0 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [A B C D] = GroupPoints(numiter, noise_var)");
		return;
	}
	SimplePoint::_id = 0;

	vector<SimpleLinkedPoint*> P;
	SimplePointFactory& factory = SimplePointFactory::getInstance();
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
			int gid = GetData2(P0, i, 3, dimsP[0], dimsP[1], (float)0);
			SimplePoint* p = factory.makePont(x, y, 0, fr, gid);
			SimpleLinkedPoint* lp = new SimpleLinkedPoint(p);
			P.push_back(lp);
		}
	}

	int numIter = 1;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[1], classMode);
	}
	float sigma = 100.0;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(sigma, prhs[2], classMode);
	}
	float rate = 0.25;
	vector<vector<SimpleLinkedPoint*>> frames = organizeIntoFrames(P);
	int numFrames = frames.size();
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(numFrames, prhs[3], classMode);
	}

	vector<SimpleLinkedPoint*> Q = frames[0];
	vector<SimpleLinkedPoint*> Z = frames[0]; //points that are seeking for a link
	int sigmaN = 1;
	for (int k = 1; k < numFrames; ++k)
	{
		printf("Frame: %d\n", k);
		Z = addFrame(Z, frames[k], sigma);
		Q.insert(Q.end(), frames[k].begin(), frames[k].end());
		updateSigma(frames[k], sigma, sigmaN);
	}
	vector<int> labels = clusterPoints(Q);

	if (nlhs >= 1)
	{
		const int dims[] = { Q.size(), 5 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SimplePoint* p = Q[i]->p;
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
		const int dims[] = { Q.size(), 5 };
		vector<int> F(dims[0] * dims[1], -1);
		for (int i = 0; i < dims[0]; ++i)
		{
			SimplePoint* p = Q[i]->p;
			SetData2(F, i, 0, dims[0], dims[1], Q[i]->prev == NULL ? -1: Q[i]->prev->p->id);
			SetData2(F, i, 1, dims[0], dims[1], Q[i]->p->id);
			SetData2(F, i, 0, dims[0], dims[1], Q[i]->next == NULL ? -1 : Q[i]->next->p->id);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	for (int i = 0; i < P.size(); ++i)
	{
		delete P[i];
	}
	factory.clean();
	mexUnlock();
}

