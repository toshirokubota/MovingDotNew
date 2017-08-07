#pragma once

#include <MotionPoint.h>

void
EstimateTransform(vector<MotionPoint*>& p1, const vector<MotionPoint*>& p2, int maxCandidates, float range);

float
CompatibilityTransform(MotionPoint* p1, MotionPoint* p2, int m, int n, float sgm);

float
InterframeCompatibilityTransform(MotionPoint* p1, MotionPoint* p2, int m, int n, float sgm);

void
UpdateTransformEstimate(vector<MotionPoint*>& points, float sigma, float rate);

void
UpdateLinkWeights(vector<MotionPoint*>& points, float thres, float sigma);

void
UpdateProbMeasureTmp(vector<MotionPoint*>& frame1, vector<MotionPoint*>& frame2, float thres, float sigma);

void
UpdateProbMeasure(vector<MotionPoint*>& frame1);

vector<int>
labelIntraFramePoints(vector<MotionPoint*>& points);

vector<pair<int, int>>
findIntraframeMatches(vector<MotionPoint*>& frame);

vector<pair<int, int>>
findInterframeMatches(vector<MotionPoint*>& frame1, vector<MotionPoint*>& frame2);
