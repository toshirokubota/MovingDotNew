#include <Canny.h>
#include <szmexutilitytemplate.h>
#include <szMexUtility.h>

const float M_PI = 3.14159265358979323846264338327;

vector<float>
CannyNonMaximumSuppression(const vector<float>& dmag,
	const vector<float>& dy,
	const vector<float>& dx,
	int ndim, const int* dims)
{
	vector<float> res(numberOfElements(ndim, dims), 0);
	for (int j = 1; j<dims[1] - 1; ++j)
	{
		for (int i = 1; i<dims[0] - 1; ++i)
		{
			float ux = GetData2(dx, i, j, dims[0], dims[1], 0.0f);
			float uy = GetData2(dy, i, j, dims[0], dims[1], 0.0f);
			float or = atan2(uy, ux);

			float gc = GetData2(dmag, i, j, dims[0], dims[1], 0.0f);
			float gw = GetData2(dmag, i - 1, j, dims[0], dims[1], 0.0f);
			float ge = GetData2(dmag, i + 1, j, dims[0], dims[1], 0.0f);
			float gn = GetData2(dmag, i, j - 1, dims[0], dims[1], 0.0f);
			float gs = GetData2(dmag, i, j + 1, dims[0], dims[1], 0.0f);
			float gsw = GetData2(dmag, i - 1, j + 1, dims[0], dims[1], 0.0f);
			float gse = GetData2(dmag, i + 1, j + 1, dims[0], dims[1], 0.0f);
			float gnw = GetData2(dmag, i - 1, j - 1, dims[0], dims[1], 0.0f);
			float gne = GetData2(dmag, i + 1, j - 1, dims[0], dims[1], 0.0f);
			float g0;
			if (j == 128 && i == 128)
			{
				i += 0;
			}
			if (ux*uy>0)
			{
				if (Abs(ux)<Abs(uy))
				{
					if ((g0 = Abs(uy*gc))<=Abs(ux*gse + (uy - ux)*gs) ||
						g0<=Abs(ux*gnw + (uy - ux)*gn)) //used to be <=
						continue;
				}
				else
				{
					if ((g0 = Abs(ux*gc))<=Abs(uy*gse + (ux - uy)*ge) ||
						g0<=Abs(uy*gnw + (ux - uy)*gw)) //used to be <=
						continue;
				}
			}
			else
			{
				if (Abs(ux)<Abs(uy))
				{
					if ((g0 = Abs(uy*gc))<Abs(ux*gne - (uy + ux)*gn) ||
						g0<=Abs(ux*gsw - (uy + ux)*gs)) //used to be <=
						continue;
				}
				else
				{
					if ((g0 = Abs(ux*gc))<=Abs(uy*gne - (ux + uy)*ge) ||
						g0<Abs(uy*gsw - (ux + uy)*gw)) //used to be <=
						continue;
				}
			}
			SetData2(res, i, j, dims[0], dims[1], gc);
		}
	}
	return res;
}

float
ThresholdSelection(const vector<float> grd,
	float percent,
	int ndim, const int* dims)
{
	int i;
	int nump = dims[0] * dims[1];
	float maxval = grd[0];
	float minval = maxval;
	for (i = 1; i<nump; ++i)
	{
		float val = grd[i];
		if (maxval<val)
			maxval = val;
		if (minval>val)
			minval = val;
	}
	if (maxval == minval) //constant image
		return maxval;

	int numbin = 256;
	vector<float> vhist(numbin, 0);
	for (i = 0; i<nump; ++i)
	{
		float val = grd[i];
		int id = (int)((numbin - 1)*(val - minval) / (maxval - minval));
		vhist[id]++;
	}
	bool found = false;
	int idx = numbin;
	do
	{
		int target = Round(nump*percent);
		int sum = 0;
		for (i = 0; i<numbin; ++i)
		{
			sum += vhist[i];
			if (sum>target)
				break;
		}
		//for image with small number of edges, (in particular synthetic one)
		//the selection may still be the same with minval.
		//if so, increase the percentile and try again.
		if (i>0)
		{
			found = true;
			idx = i;
		}
		else
			percent = 1 - (1 - percent) / 2;
	} while (!found);
	return (float)idx*(maxval - minval) / numbin + minval;
}

void
follow(vector<float>& trace,
	const vector<float>& edge,
	vector<unsigned char>& flag,
	int y, int x,
	float low,
	const int* dims)
{
	SetData2(flag, x, y, dims[0], dims[1], (unsigned char)1);
	SetData2(trace, x, y, dims[0], dims[1], GetData2(edge, x, y, dims[0], dims[1], 0.0f));
	for (int i = y - 1; i <= y + 1; ++i)
	{
		for (int j = x - 1; j <= x + 1; ++j)
		{
			if (!GetData2(flag, j, i, dims[0], dims[1], (unsigned char)1))
			{
				if (GetData2(edge, j, i, dims[0], dims[1], (float)0)>low)
				{
					follow(trace, edge, flag, i, j, low, dims);
				}
			}
		}
	}
}

vector<float>
EdgeTrace(const vector<float>& edge,
	float high, float low,
	const int* dims)
{
	vector<float> res(dims[0] * dims[1], .0f);
	vector<unsigned char> flag(dims[0] * dims[1], (unsigned char)0);
	for (int i = 0; i<dims[1]; ++i)
	{
		for (int j = 0; j<dims[0]; ++j)
		{
			if (!GetData2(flag, j, i, dims[0], dims[1], (unsigned char)1) && GetData2(edge, j, i, dims[0], dims[1], (float)0) > high)
			{
				follow(res, edge, flag, i, j, low, dims);
			}
		}
	}
	return res;
}

#include <szParticle.h>
#include <algorithm>
#include <set>
using namespace std;

/*
From P, trace edges in the direction initially set in (DIRX, DIRY).
The trace continues until the edge strength goes down below LOW or
it meets an edge that has already been processed.
*/
vector<CParticle*> 
_track(CParticle* p, 
	vector<CParticle*>& mp,
	vector<bool>& sm,
	float dirX, float dirY, float low,
	int ndim, const int* dims)
{
	vector<CParticle*> tr;
	while (p != NULL)
	{
		if (p->m_Life < low) break;
		if (GetData2(sm, p->m_X, p->m_Y, dims[0], dims[1], true)) break;
		SetData2(sm, p->m_X, p->m_Y, dims[0], dims[1], true);
		tr.push_back(p);
		CParticle* q = NULL;
		float dotmax = 0;
		for (int i = -1; i <= 1; ++i)
		{
			for (int j = -1; j <= 1; ++j)
			{
				if (i == 0 && j == 0) continue;
				CParticle* r = GetData2(mp, p->m_X + j, p->m_Y + i, dims[0], dims[1], (CParticle*)NULL);
				if (r != NULL)
				{
					float dotval = j*dirX + i*dirY;
					if (dotval > dotmax)
					{
						dotmax = dotval;
						q = r;
					}
				}
			}
		}
		if (q != NULL)
		{
			dirX = q->m_X - p->m_X;
			dirY = q->m_Y - p->m_Y;
			p = q;
		}
	}
	return tr;
}
#include <mex.h>

vector<float>
EdgeTrace(const vector<float>& mag,
	const vector<float>& gy,
	const vector<float>& gx,
	float high, float low,
	int ndim, const int* dims)
{
	vector<pair<float,CParticle*>> pnts;
	vector<CParticle*> mp(numberOfElements(ndim, dims), NULL);
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			float val = GetData2(mag, j, i, dims[0], dims[1], 0.0f);
			if (val > 0.0f)
			{
				CParticle* p = new CParticle(j, i, 0, val);
				SetData2(mp, j, i, dims[0], dims[1], p);
				if (val >= high)
				{
					pnts.push_back(pair<float, CParticle*>(val, p));
				}
			}
		}
	}
	sort(pnts.begin(), pnts.end());
	vector<bool> state(mp.size(), false);
	vector<float> edge(mp.size(), 0);
	for(int idx=pnts.size()-1; idx>=0; idx--)
	{
		CParticle* p = pnts[idx].second;
		if (GetData2(state, p->m_X, p->m_Y, dims[0], dims[1], true)) continue;

		vector<CParticle*> track;
		float dx = GetData2(gx, p->m_X, p->m_Y, dims[0], dims[1], 0.0f);
		float dy = GetData2(gy, p->m_X, p->m_Y, dims[0], dims[1], 0.0f);
		float dir1 = atan2(-dx, dy); //along the gradient
		float cs1 = cos(dir1), sn1 = sin(dir1);
		vector<CParticle*> tr = _track(p, mp, state, cs1, sn1, low, ndim, dims);

		//reset the state of p so that we can continue on the other direction
		SetData2(state, p->m_X, p->m_Y, dims[0], dims[1], false);
		float dir2 = atan2(dx, dy);  //the opposite direction
		float cs2 = cos(dir2), sn2 = sin(dir2);
		vector<CParticle*> tr2 = _track(p, mp, state, cs2, sn2, low, ndim, dims);
		//append tr2 in the reverse order, excluding the first one (i.e. p)
		for (int i = 1; i < tr2.size(); ++i)
		{
			tr.insert(tr.begin(), tr2[i]);
		}
		for (int i = 0; i < tr.size(); ++i)
		{
			SetData2(state, tr[i]->m_X, tr[i]->m_Y, dims[0], dims[1], true);
		}
		if (tr.size() >= 10)
		{
			for (int i = 0; i < tr.size(); ++i)
			{
				//SetData2(edge, tr[i]->m_X, tr[i]->m_Y, dims[0], dims[1], p->m_Life);
				SetData2(edge, tr[i]->m_X, tr[i]->m_Y, dims[0], dims[1], (float)(pnts.size() - idx));
			}
			printf("Edge %d: %f %d\n", pnts.size() - idx, p->m_Life, tr.size());
		}
	}
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i]) delete mp[i];
	}

	return edge;
}


/*
Morphological thinning based on
Gua & Hall, CACM, 1989, p359-373.
One operation is enough after the non-maximum suppression.
*/
vector<unsigned char>
CannyThinning(const vector<unsigned char>& image,
	const int* dims)
{
	const int NumNeighbors = 8;
	int ax[NumNeighbors] = { 1,1,0,-1,-1,-1,0,1 };
	int ay[NumNeighbors] = { 0,-1,-1,-1,0,1,1,1 };

	vector<unsigned char> res = image;
	int i, j, m;
	for (int pass = 0; pass<2; ++pass)
	{
		for (i = 1; i<dims[1] - 1; ++i)
		{
			for (j = 1; j<dims[0] - 1; ++j)
			{
				if (GetData2(res, j, i, dims[0], dims[1], (unsigned char)0))
				{
					int cp = 0;
					int np1 = 0, np2 = 0;
					for (m = 1; m <= NumNeighbors / 2; m++) {
						int y1 = i + ay[(2 * m - 2) % NumNeighbors];
						int x1 = j + ax[(2 * m - 2) % NumNeighbors];
						int y2 = i + ay[(2 * m - 1) % NumNeighbors];
						int x2 = j + ax[(2 * m - 1) % NumNeighbors];
						int y3 = i + ay[(2 * m) % NumNeighbors];
						int x3 = j + ax[(2 * m) % NumNeighbors];
						unsigned char p1 = GetData2(res, x1, y1, dims[0], dims[1], (unsigned char)0);
						unsigned char p2 = GetData2(res, x2, y2, dims[0], dims[1], (unsigned char)0);
						unsigned char p3 = GetData2(res, x3, y3, dims[0], dims[1], (unsigned char)0);
						if (!p1 && (p2 || p3))
							cp++;
						if (p1 || p2)
							np1++;
						if (p2 || p3)
							np2++;
					}
					int np = Min(np1, np2);
					if (cp == 1 && np >= 2 && np <= 3)
					{
						if (pass == 0)
						{
							unsigned char p2 = GetData2(res, j + ax[1], i + ay[1], dims[0], dims[1], (unsigned char)0);
							unsigned char p3 = GetData2(res, j + ax[2], i + ay[2], dims[0], dims[1], (unsigned char)0);
							unsigned char p8 = GetData2(res, j + ax[7], i + ay[7], dims[0], dims[1], (unsigned char)0);
							unsigned char p1 = GetData2(res, j + ax[0], i + ay[0], dims[0], dims[1], (unsigned char)0);
							if (!((p2 || p3 || !p8) && p1))
							{
								SetData2(res, j, i, dims[0], dims[1], (unsigned char)0);
							}
						}
						else if (pass == 1)
						{
							unsigned char p4 = GetData2(res, j + ax[3], i + ay[3], dims[0], dims[1], (unsigned char)0);
							unsigned char p5 = GetData2(res, j + ax[4], i + ay[4], dims[0], dims[1], (unsigned char)0);
							unsigned char p6 = GetData2(res, j + ax[5], i + ay[5], dims[0], dims[1], (unsigned char)0);
							unsigned char p7 = GetData2(res, j + ax[6], i + ay[6], dims[0], dims[1], (unsigned char)0);
							if (!((p6 || p7 || !p4) && p5))
							{
								SetData2(res, j, i, dims[0], dims[1], (unsigned char)0);
							}
						}
					}
				}
			}
		}
	}
	return res;
}

//#include <mex.h>

#include <szConvolutionTemplate.h>
void
CannyEdge(vector<unsigned char>& edge,
		vector<float>& im,
		float lowThres, float highThres, float sgm,
		int ndim, const int* dims)
{
	//vector<float> filter = makeGaussianFilter(0.0f, sgm, (int)(2 * sgm) * 2 + 1, true);
	vector<float> sm(numberOfElements(ndim, dims), 0.0f);
	//separableConvolution(sm, im, filter, CBE_MIRROR, ndim, dims);
	sm = im;
	vector<float> gr(numberOfElements(ndim, dims), 0.0f);
	vector<float> dx(numberOfElements(ndim, dims), 0.0f);
	vector<float> dy(numberOfElements(ndim, dims), 0.0f);
	for (int i = 1; i < dims[1] - 1; ++i)
	{
		for (int j = 1; j < dims[0] - 1; ++j)
		{
			float nval = GetData2(sm, j, i - 1, dims[0], dims[1], 0.0f);
			float sval = GetData2(sm, j, i + 1, dims[0], dims[1], 0.0f);
			float wval = GetData2(sm, j - 1, i, dims[0], dims[1], 0.0f);
			float eval = GetData2(sm, j + 1, i, dims[0], dims[1], 0.0f);
			float dfx = (eval - wval) / 2;
			float dfy = (sval - nval) / 2;
			SetData2(dx, j, i, dims[0], dims[1], dfx);
			SetData2(dy, j, i, dims[0], dims[1], dfy);
			SetData2(gr, j, i, dims[0], dims[1], sqrt(dfx*dfx + dfy*dfy));
		}
	}
	float percent = 0.7f;
	float ratio = 0.4f;
	vector<float> emag = CannyNonMaximumSuppression(gr, dy, dx, ndim, dims);
	emag = EdgeTrace(emag, highThres, lowThres, dims);

	for (int i = 0; i<dims[0] * dims[1]; ++i)
	{
		edge[i] = emag[i] > 0 ? 255 : 0;
	}

	edge = CannyThinning(edge, dims);
	edge = CannyThinning(edge, dims);
}

void
CannyEdge(vector<unsigned char>& edge,
	vector<float>& im, //input - no Gaussian smoothing included.
	int ndim, const int* dims)
{
	vector<float> gr(numberOfElements(ndim, dims), 0.0f);
	vector<float> dx(numberOfElements(ndim, dims), 0.0f);
	vector<float> dy(numberOfElements(ndim, dims), 0.0f);
	for (int i = 1; i < dims[1] - 1; ++i)
	{
		for (int j = 1; j < dims[0] - 1; ++j)
		{
			float val = GetData2(im, j, i, dims[0], dims[1], 0.0f);
			float nval = GetData2(im, j, i - 1, dims[0], dims[1], val);
			float sval = GetData2(im, j, i + 1, dims[0], dims[1], val);
			float wval = GetData2(im, j - 1, i, dims[0], dims[1], val);
			float eval = GetData2(im, j + 1, i, dims[0], dims[1], val);
			float dfx = (eval - wval) / 2;
			float dfy = (sval - nval) / 2;
			SetData2(dx, j, i, dims[0], dims[1], dfx);
			SetData2(dy, j, i, dims[0], dims[1], dfy);
			SetData2(gr, j, i, dims[0], dims[1], sqrt(dfx*dfx + dfy*dfy));
		}
	}
	float percent = 0.7f;
	float ratio = 0.4f;
	vector<float> emag = CannyNonMaximumSuppression(gr, dy, dx, ndim, dims);
	float high = ThresholdSelection(gr, percent, ndim, dims);
	float low = high*ratio;
	printf("low thres = %f, high thres = %f\n", low, high);
	//emag = EdgeTrace(emag, high, low, dims);
	emag = EdgeTrace(emag, dy, dx, high, low, ndim, dims);

	//vector<unsigned char> edge(dims[0] * dims[1], (unsigned char)0);
	for (int i = 0; i<dims[0] * dims[1]; ++i)
	{
		edge[i] = emag[i] > 0 ? 255 : 0;
	}

	edge = CannyThinning(edge, dims);
	edge = CannyThinning(edge, dims);
}
