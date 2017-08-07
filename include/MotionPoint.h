#ifndef _MotionPoint_h_
#define _MotionPoint_h_

#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;
#include <szMexUtility.h>

struct MotionPoint
{
	MotionPoint(float x0 = 0, float y0 = 0)
	{
		x = x0;
		y = y0;
		id = _id++;
	}

	void setCandidates(vector<MotionPoint*>& cand)
	{
		candidates = cand;
		for (int i = 0; i < cand.size(); ++i)
		{
			MotionPoint* c = cand[i];
			velocityX.push_back(c->x - x);
			velocityY.push_back(c->y - y);
			probCor.push_back(1.0 / (cand.size() + 1));
		}
		probCor.push_back(1.0 / (cand.size() + 1));
		probCor0 = probCor;
		velocityX0 = velocityX;
		velocityY0 = velocityY;
	}

	void initializeLinks(vector<MotionPoint*>& frame)
	{
		vector<vector<float>> zeros();
		probLink = vector<vector<vector<float>>>(candidates.size());
		probLink0 = vector<vector<vector<float>>>(candidates.size());
		for (int c = 0; c < candidates.size(); ++c)
		{
			probLink[c] = vector<vector<float>>(frame.size());
			probLink0[c] = vector<vector<float>>(frame.size());
			for (int i = 0; i < frame.size(); ++i)
			{
				probLink[c][i] = vector<float>(frame[i]->candidates.size(), 0.5f);
				probLink0[c][i] = vector<float>(frame[i]->candidates.size(), 0.5f);
			}
		}
	}
	float getLink(int cand, int pixel, int cand2)
	{
		return probLink[cand][pixel][cand2];
	}
	void setLink0(int cand, int pixel, int cand2, float val)
	{
		probLink0[cand][pixel][cand2] = val;
	}
	void updateLinks()
	{
		for (int i = 0; i < probLink.size(); ++i)
		{
			for (int j = 0; j < probLink[i].size(); ++j)
			{
				for (int k = 0; k < probLink[i][j].size(); ++k)
				{
					probLink[i][j][k] = probLink0[i][j][k];
				}
			}
		}
	}

	vector<float> probCor; //inter-framce correspondence probability
	vector<vector<vector<float>>> probLink; //intra-frame grouping probability - AWFUL AWFUL 3D vectors...
	vector<float> velocityX; //velocity estimate
	vector<float> velocityY; //velocity estimate

	vector<MotionPoint*> candidates; //possible candidates for inter-frame correspondence
	float x;
	float y;

	vector<float> probCor0; //temporary copy of inter-framce correspondence probability
	vector<vector<vector<float>>> probLink0; //temporary copy of intra-frame grouping probability - AWFUL AWFUL 3D vectors...
	vector<float> velocityX0; //temporary copy of velocity estimate
	vector<float> velocityY0; //temporary copy of velocity estimate

	int id;
	static int _id;
};

struct MotionPointSynth : public MotionPoint
{
	MotionPointSynth(float x0 = 0, float y0 = 0, float ux=0, float uy=0, int gid=0) : MotionPoint(x0, y0)
	{
		this->ux = ux;
		this->uy = uy;
		groupId = gid;
	}
	void Move(float rate, float var)
	{
		y += rate*(uy + var*(rndm(0) - .5));
		x += rate*(ux + var*(rndm(0) - .5));
	}
	float ux;
	float uy;
	int groupId;
};

#endif /* _MotionPoint_H */
