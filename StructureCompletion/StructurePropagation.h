#pragma  once
#include "OpenCvUtility.h"
#include "PointManager.h"
#include <map>
#include <list>
#include <set>
#include <math.h>

class StructurePropagation
{
public:
	StructurePropagation() {
		ki = 0.5;
		ks = 0.5;
	}
	~StructurePropagation(){}
	void Run(const Mat1b &_mask, const Mat& _img, const vector<vector<Point>> &linePoints, Mat& result);
	//利用目前类中已经存储的数据继续经行修补
	//void RunAgain();
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
	double calcEs(PointPos i, PointPos xi);
	double calcEi(const Mat &mat, PointPos i, PointPos xi);
	double calcE2(const Mat &mat, PointPos i1, PointPos i2);
	void calcMij(Node &n, const Mat &mat, const vector<PointPos> &samplePoints);
};