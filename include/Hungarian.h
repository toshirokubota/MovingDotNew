#pragma once
#include <vector>
using namespace std;

struct Hungarian {
	/* 
	C is a cost matrix as a flat vector, and dims is the dimension of the matrix.
	The first dimension is considered as jobs, and the second dimension is considered as Workers.
	*/
	Hungarian(vector<float>& C, const int* dims)
	{
		this->C = C;
		nJobs = dims[0];
		nWorkers = dims[1];
		jobs = vector<int>(nWorkers, -1); //map worker to job
		workers = vector<int>(nJobs, -1); //map job to worker
	}
	float solve(bool bMaximize = true);
	int getWorker(int i) const
	{
		return workers[i];
	}	
	int getJob(int i) const
	{
		return jobs[i];
	}
	float getCost(int i, int j) const
	{
		return C[i*nJobs + j];
	}
	int nJobs;
	int nWorkers;
	vector<float> C;
	vector<int> jobs; //map job to worker after matching.
	vector<int> workers; //map worker to job after matching.
};