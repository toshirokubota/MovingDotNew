#include <MotionPoint.h>
#include <DisjointSet.h>

float
_Dist(float y1, float x1, float y2, float x2)
{
	return (y1 - y2)*(y1 - y2) + (x1 - x2)*(x1 - x2);
}

void
EstimateTransform(vector<MotionPoint*>& p1, const vector<MotionPoint*>& p2, int maxCandidates, float range)
{
	int nump = p1.size();
	for (int i = 0; i<nump; ++i)
	{
		float x0 = p1[i]->x;
		float y0 = p1[i]->y;
		int count = 0;
		vector<pair<float, MotionPoint*>> pairs;
		for (int j = 0; j<nump; ++j)
		{
			float x1 = p2[j]->x;
			float y1 = p2[j]->y;
			float dd = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
			pairs.push_back(pair<float, MotionPoint*>(dd, p2[j]));
		}
		sort(pairs.begin(), pairs.end());
		vector<MotionPoint*> candidates;
		for (int j = 0; j < maxCandidates; ++j)
		{
			if (pairs[j].first > 2 * range * range)
			{
				break;
			}
			candidates.push_back(pairs[j].second);
		}
		p1[i]->setCandidates(candidates);
	}
}

float
CompatibilityTransform(MotionPoint* p1, MotionPoint* p2, int m, int n, float sgm)
{
	float vy1 = p1->velocityY[m];
	float vx1 = p1->velocityX[m];
	float vy2 = p2->velocityY[n];
	float vx2 = p2->velocityX[n];
	float dd = _Dist(vy1, vx1, vy2, vx2);
	if (dd > 2.*sgm*sgm)
	{
		return .0;
	}
	else
	{
		return exp(-dd / (2.*sgm*sgm));
	}
}

float
InterframeCompatibilityTransform(MotionPoint* p1, MotionPoint* p2, int m, int n, float sgm)
{
	float vy1 = p1->velocityY[m];
	float vx1 = p1->velocityX[m];
	float vy2 = p2->velocityY[n];
	float vx2 = p2->velocityX[n];
	float dd = _Dist(vy1, vx1, vy2, vx2);
	if (dd > 2.*sgm*sgm)
	{
		return .0;
	}
	else
	{
		return exp(-dd / (2.*sgm*sgm));
	}
}

void
UpdateTransformEstimate(vector<MotionPoint*>& points, float sigma, float rate)
{
	int nump = points.size();
	int i, m, n;
	float alpha = 0.5;
	for (i = 0; i < nump; ++i)
	{
		float total_sum = 0;
		MotionPoint* pnt1 = points[i];
		int numc = pnt1->candidates.size();
		vector<float> sum(numc, 0.0f);
		for (m = 0; m < pnt1->candidates.size(); ++m)
		{
			float evy1 = pnt1->velocityY[m];
			float evx1 = pnt1->velocityX[m];
			float nvx = 0.0f;
			float nvy = 0.0f;
			for (int j = 0; j < nump; ++j)
			{
				if (i != j)
				{
					MotionPoint* pnt2 = points[j];
					for (n = 0; n < pnt2->candidates.size(); ++n)
					{
						float evy2 = pnt2->velocityY[n];
						float evx2 = pnt2->velocityX[n];
						float lnk1 = pnt1->getLink(m, j, n); // probLink[j][m][n];
						float lnk2 = pnt2->getLink(n, i, m); // probLink[i][n][m];
						float dd = CompatibilityTransform(pnt1, pnt2, m, n, sigma);
						float dy = evy2 - evy1;
						float dx = evx2 - evx1;
						float pp = pnt2->probCor[n];
						nvy += dd*pp*lnk1*lnk2*dy;
						nvx += dd*pp*lnk1*lnk2*dx;
					}
				}
			}

			pnt1->velocityY0[m] = evy1 + rate*nvy;
			pnt1->velocityX0[m] = evx1 + rate*nvx;
		}
	}
	//copy back from the temporary velocity
	for (i = 0; i < nump; ++i)
	{
		MotionPoint* pnt1 = points[i];
		for (m = 0; m < pnt1->candidates.size(); ++m)
		{
			pnt1->velocityX[m] = pnt1->velocityX0[m];
			pnt1->velocityY[m] = pnt1->velocityY0[m];
		}
	}
}

void
UpdateLinkWeights(vector<MotionPoint*>& points, float thres, float sigma)
{
	int nump = points.size();
	for (int i = 0; i < nump; ++i)
	{
		MotionPoint* p1 = points[i];
		int numc1 = p1->candidates.size();
		for (int m = 0; m < numc1; ++m)
		{
			for (int j = 0; j < nump; ++j)
			{
				MotionPoint* p2 = points[j];
				int numc2 = p2->candidates.size();
				for (int n = 0; n < numc2; ++n)
				{
					if (i != j)
					{
						float dd = CompatibilityTransform(p1, p2, m, n, sigma);
						float lnk1 = p1->getLink(m, j, n); // probLink[j][m][n];
						float lnk2 = p2->getLink(n, i, m); // [i][n][m];
						float pp = p2->probCor[n];
						float ee = pp*lnk1*lnk2;
						float pp2 = ee*dd / (ee*dd + (1. - ee)*thres);
						if (pp2 == pp2)
						{
							//res[i].SetLink(j, m, n, pp2);
							p1->setLink0(m, j, n, pp2); // probLink0[j][m][n] = pp2;
						}
						else
						{
							//res[i].SetLink(j, m, n, 0.5f);
							p1->setLink0(m, j, n, 0.5f); // probLink0[j][m][n] = 0.5f;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < nump; ++i)
	{
		points[i]->updateLinks();
	}
}

void
UpdateProbMeasureTmp(vector<MotionPoint*>& frame1, vector<MotionPoint*>& frame2,
float thres, float sigma)
{
	int nump = frame1.size();
	for (int i = 0; i < nump; ++i)
	{
		MotionPoint* p1 = frame1[i];
		int numc1 = p1->candidates.size();
		float total_sum = 0;
		vector<float> sum(numc1, 0.0f);
		for (int m = 0; m < numc1; ++m)
		{ //compute fitness for each candidate
			MotionPoint* p3 = p1->candidates[m];
			float match_prob = .0;
			//check if the candidate also consider it as a candidate
			for (int m2 = 0; m2 < p3->candidates.size(); ++m2)
			{
				if (p3->candidates[m2] == p1)
				{
					match_prob = p3->probCor[m2];
					break;
				}
			}
			//if(match_prob>0) 
			{ //mutual candidancy
				for (int j = 0; j < nump; ++j)
				{
					//compute the support from other INTRA-frame points
					//measured by 1. transformation similarity
					//            2. strength of links
					if (i != j)
					{
						MotionPoint* p2 = frame1[j];
						int numc2 = p2->candidates.size();
						for (int n = 0; n < numc2; ++n)
						{
							float dd = CompatibilityTransform(p1, p2, m, n, sigma);
							float lnk1 = p1->getLink(m, j, n); // probLink[j][m][n]; // p1.GetLink(j, m, n)
							float lnk2 = p2->getLink(n, i, m); // probLink[i][n][m]; // p2.GetLink(i, n, m);
							float pp = p2->probCor[n]; // p2.GetProb(n);
							float ss = pp*lnk1*lnk2;
							sum[m] += ss*match_prob;
						}
					}
				}
				sum[m] *= p1->probCor[m]; // p1.GetProb(m);
				total_sum += sum[m];
			}
		}
		if (total_sum > 1.0e-5)
		{
			for (int m = 0; m < numc1; ++m)
			{
				//res[i].SetProb(m, sum[m] / (total_sum + p1.GetProb(numc1)*thres));
				p1->probCor0[m] = sum[m] / (total_sum + p1->probCor[numc1] * thres);
			}
			//res[i].SetProb(numc1, p1.GetProb(numc1)*thres / (total_sum + p1.GetProb(numc1)*thres));
			p1->probCor0[numc1] = p1->probCor[numc1] * thres / (total_sum + p1->probCor[numc1] * thres);
		}
		else
		{
			for (int m = 0; m < numc1; ++m)
			{
				//res[i].SetProb(m, 0.0f);
				p1->probCor0[m] = 0.0f;
			}
			//res[i].SetProb(numc1, 1.0f);
			p1->probCor0[numc1] = 1.0f;
		}
	}
}

void
UpdateProbMeasure(vector<MotionPoint*>& frame1)
{
	int nump = frame1.size();
	for (int i = 0; i < nump; ++i)
	{
		MotionPoint* p1 = frame1[i];
		for (int m = 0; m < p1->probCor.size(); ++m)
		{
			p1->probCor[m] = p1->probCor0[m];
		}
	}
}

vector<int>
labelIntraFramePoints(vector<MotionPoint*>& points)
{
	vector<Node<int>*> vnodes;
	for (int i = 0; i < points.size(); ++i)
	{
		vnodes.push_back(makeset(i));
	}
	for (int i = 0; i < points.size(); ++i)
	{
		for (int j = i + 1; j < points.size(); ++j)
		{
			//if (i == j) continue;
			for (int m = 0; m < points[i]->candidates.size(); ++m)
			{
				for (int n = 0; n<points[j]->candidates.size(); ++n)
				{
					//if (points[i].GetLink(j, m, n) > .5) {
					if (points[i]->getLink(m, j, n) > 0.5) //probLink[j][m][n] > .5)
					{
						merge(vnodes[i], vnodes[j]);
					}
				}
			}
		}
	}
	vector<Node<int>*> labels = clusters(vnodes);
	vector<int> ilabels(points.size());
	for (int i = 0; i < points.size(); ++i)
	{
		ilabels[i] = distance(labels.begin(), find(labels.begin(), labels.end(), findset(vnodes[i])));
	}
	for (int i = 0; i < points.size(); ++i)
	{
		delete vnodes[i];
	}
	return ilabels;
}

vector<pair<int, int>>
findIntraframeMatches(vector<MotionPoint*>& frame)
{
	vector<pair<int, int>> matches;
	for (int i = 0; i < frame.size(); ++i)
	{
		for (int j = 0; j < frame.size(); ++j)
		{
			for (int m = 0; m < frame[i]->candidates.size(); ++m)
			{
				for (int n = 0; n < frame[j]->candidates.size(); ++n)
				{
					//float w = frame[i].GetLink(j, m, n);
					float w = frame[i]->getLink(m, j, n); // probLink[j][m][n];
					if (w > 0.5)
					{
						matches.push_back(pair<int, int>(i, j));
						break;
					}
				}
			}
		}
	}
	return matches;
}

vector<pair<int, int>>
findInterframeMatches(vector<MotionPoint*>& frame1, vector<MotionPoint*>& frame2)
{
	vector<pair<int, int>> matches;
	for (int i = 0; i < frame1.size(); ++i)
	{
		for (int m = 0; m < frame1[i]->candidates.size(); ++m)
		{
			float p = frame1[i]->probCor[m];
			if (p > 0.5)
			{
				//int j = frame1[i].GetIndex(m);
				MotionPoint* p2 = frame1[i]->candidates[m];
				for (int n = 0; n < p2->candidates.size(); ++n)
				{
					float q = p2->probCor[n];
					if (q > 0.5 && p2->candidates[n] == frame1[i])
					{
						matches.push_back(pair<int, int>(frame1[i]->id, p2->id));
						break;
					}
				}
			}
		}
	}
	return matches;
}
