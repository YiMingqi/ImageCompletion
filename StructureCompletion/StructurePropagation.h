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
		ki = 0.75;
		ks = 0.25;
		//ki = 0.25;
		//ks = 0.75;
	}
	~StructurePropagation(){}
	void Run(const Mat1b &_mask, const Mat& _img, Mat1b &Linemask, vector<vector<Point>> &linePoints, Mat& result);
	void StructurePropagation::TextureCompletion(const Mat1b &_mask, Mat1b &LineMask, const Mat &mat, Mat &result);
	void StructurePropagation::TextureCompletion2(Mat1b _mask, Mat1b LineMask, const Mat &mat, Mat &result);
	void SetParm(int _blocksize,int _samplestep,int _iscurve);

private:
	int blockSize;
	int sampleStep;
	int isCurve;
	double ki;
	double ks;
	PointManager pointManager;

	void getResult(Mat1b mask, int *sampleIndices, const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, Mat& result);
	void ModifyMask(Mat1b &LineMask, vector<PointPos>AnchorPoints);
	int *DP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat);
	int *BP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat);
	double gauss(double x);
	double calcEs(const PointPos &i, const PointPos &xi);
	double calcEi(const Mat &mat, const PointPos &i, const PointPos &xi);
	double calcE2(const Mat &mat, const PointPos &i1, const PointPos &i2, const PointPos &xi1, const PointPos &xi2);
	void calcMij(Node &n, const list<shared_ptr<Edge>>::iterator &edgeItor, const Mat &mat, const vector<PointPos> &samplePoints);
};