#include <MovingEdgeDetector.h>
#include <szMexUtility.h>
#include <szmexutilitytemplate.h>

namespace TK
{
	void
	MovingEdgeDetector::run(vector<float>& im, int ndim, const int* dims)
	{
		//reset the edges
		for (int i = 0; i < mp.size(); ++i)
		{
			if (mp[i]) delete mp[i];
		}
		mp = vector<CParticle*>(numberOfElements(ndim, dims), NULL);
		edges.clear();

		this->ndim = ndim;
		this->dims[0] = dims[0];
		this->dims[1] = dims[1];

		Canny* prev = latestDetection;
		edgeDetector[next]->run(im, ndim, dims);
		latestDetection = edgeDetector[next];
		valid[next] = true;

		for (int i = 0; i < latestDetection->edges.size(); ++i)
		{
			int count = 0;
			int numEdges = latestDetection->edges[i].size();
			for (int j = 0; j < numEdges; ++j)
			{
				CParticle* p = latestDetection->edges[i][j];
				for (int n = 0; n < NumDetectors; ++n)
				{
					if (latestDetection == edgeDetector[n]) continue;
					if (!valid[n]) continue;

					if (GetData2(edgeDetector[n]->mp, p->m_X, p->m_Y, dims[0], dims[1], (CParticle*)NULL) == NULL)
					{
						count++;
					}
				}
			}
			if ((float)count / (float)(numEdges * (NumDetectors - 1)) > threshold)
			{
				vector<CParticle*> tr;
				for (int j = 0; j < edgeDetector[next]->edges[i].size(); ++j)
				{
					CParticle* p = edgeDetector[next]->edges[i][j];
					CParticle* q = new CParticle(*p);
					tr.push_back(q);
				}
				edges.push_back(tr);
			}
		}
		next = (next + 1) % NumDetectors;
	}
	vector<float>
	MovingEdgeDetector::retrieveMovingEdgeStrength()
	{
		vector<float> em(mp.size(), 0.0f);
		for (int i = 0; i < edges.size(); ++i)
		{
			for (int j = 0; j < edges[i].size(); ++j)
			{
				CParticle* p = edges[i][j];
				SetData2(em, p->m_X, p->m_Y, dims[0], dims[1], p->m_Life);
			}
		}
		return em;
	}
	vector<float>
	MovingEdgeDetector::retrieveEdgeStrength()
	{
		vector<float> ed(mp.size(), 0.0f);
		int index = (next - 1 + NumDetectors) % NumDetectors;
		if (valid[index])
		{
			ed = edgeDetector[index]->retrieve();
		}
		return ed;
	}
	vector<unsigned char>
	MovingEdgeDetector::retrieveMovingEdges()
	{
		vector<unsigned char> em(mp.size(), 0);
		for (int i = 0; i < edges.size(); ++i)
		{
			for (int j = 0; j < edges[i].size(); ++j)
			{
				CParticle* p = edges[i][j];
				SetData2(em, p->m_X, p->m_Y, dims[0], dims[1], (unsigned char)255);
			}
		}
		return em;
	}
	vector<unsigned char>
	MovingEdgeDetector::retrieveEdges()
	{
		vector<unsigned char> ed(mp.size(), 0.0f);
		int index = (next - 1 + NumDetectors) % NumDetectors;
		if (valid[index])
		{
			vector<float> em = edgeDetector[index]->retrieve();
			for (int i = 0; i < em.size(); ++i)
			{
				if (em[i] > 0) ed[i] = 255;
			}
		}
		return ed;
	}
}
