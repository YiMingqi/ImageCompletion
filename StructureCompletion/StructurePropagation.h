#pragma  once
#include "OpenCvUtility.h"
#include "PointManager.h"
#include "OpenCvUtility.h"
#include <map>
#include <list>
#include <set>
#include <math.h>

class StructurePropagation
{
public:
	StructurePropagation() {
		ki = 0.7;
		ks = 0.3;
	}
	~StructurePropagation(){}
	void Run(const Mat1b &_mask, const Mat& _img, vector<vector<Point>> &linePoints, Mat& result);
	void SetParm(int _blocksize,int _samplestep,int _iscurve);

private:
	int blockSize;
	int sampleStep;
	int isCurve;
	double ki;
	double ks;
	PointManager pointManager;


	void getResult(int *sampleIndices, const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, Mat& result);
	int *DP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat);
	int *BP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat);
	double calcEs(const PointPos &i, const PointPos &xi);
	double calcEi(const Mat &mat, const PointPos &i, const PointPos &xi);
	double calcE2(const Mat &mat, const PointPos &i1, const PointPos &i2, const PointPos &xi1, const PointPos &xi2);
	void calcMij(Node &n, const list<shared_ptr<Edge>>::iterator &edgeItor, const Mat &mat, const vector<PointPos> &samplePoints);
};