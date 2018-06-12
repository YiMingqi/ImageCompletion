#include "stdafx.h"
#include "StructurePropagation.h"


void StructurePropagation::Run(const Mat1b &_mask, const Mat& _img, const vector<vector<Point>> &linePoints, Mat& result) {

	// Calculate gray map
	Mat grayMat = Mat::zeros(_img.rows, _img.cols, CV_8UC1);
	for (int i = 0; i < _img.rows; i++) {
		for (int j = 0; j < _img.cols; j++) {
			Vec3b tmp = _img.at<Vec3b>(i, j);
			grayMat.at<uchar>(i, j) = (uchar)((114 * tmp[0] + 587 * tmp[1] + 299 * tmp[2] + 500) / 1000);
		}
	}

	pointManager.reset(linePoints, grayMat, blockSize);


	curLine = 0;

	//getSamples();

	/*curLine = 0;
	int *sampleIndices = NULL;
	if (lines.size() == 1) {
		sampleIndices = DP(grayMat, points);
		_img.copyTo(result);
	}
	else if (lines.size() > 1) {
		// multiple line or curve
	}*/
}

void StructurePropagation::getResult(int *sampleIndices, Mat& result) {
	int offset1 = blockSize / 2;
	int offset2 = blockSize - offset1;
	for (int i = 0; i < anchorPoints.size(); i++) {
		Point src = pointManager.getPoint(samplePoints[sampleIndices[i]]);
		Point tar = pointManager.getPoint(anchorPoints[i]);
		for (int m = -offset1; m < offset2; m++) {
			int tary = tar.y + m;
			const Vec3b* srcPtr = result.ptr<Vec3b>(src.y + m);
			for (int n = -offset1; n < offset2; n++) {
				result.at<Vec3b>(tar.y + m, tar.x + n) = srcPtr[src.x + n];
			}
		}
	}
	free(sampleIndices);
}

int *StructurePropagation::DP(const Mat &mat) {
	double *M = (double *)malloc(2 * samplePoints.size() * sizeof(double));
	int *record = (int *)malloc(samplePoints.size() * anchorPoints.size() * sizeof(int));

	

	for (int i = 0; i < samplePoints.size(); i++) {
		M[i] = ks * calcEs(anchorPoints[0], samplePoints[i]) +
			ki * calcEi(mat, anchorPoints[0], samplePoints[i]);
	}

	int i, curOffset, preOffset;
	for (i = 1; i < anchorPoints.size(); i++) {
		curOffset = (i % 2) * samplePoints.size();
		preOffset = ((i + 1) % 2) * samplePoints.size();
		for (int j = 0; j < samplePoints.size(); j++) {
			double E1 = ks * calcEs(anchorPoints[i], samplePoints[j]) +
				ki * calcEi(mat, anchorPoints[i], samplePoints[j]);
			double min = INT_MAX;
			for (int k = 0; k < samplePoints.size(); k++) {
				double tmp = calcE2(mat, samplePoints[j], samplePoints[k]) + M[preOffset + k];
				if (tmp < min) {
					record[samplePoints.size()*i + j] = k;
					min = tmp;
				}
			}
			M[curOffset + j] = E1 + min;
		}
	}

	int *sampleIndices = (int*)malloc(anchorPoints.size() * sizeof(int));
	double min = INT_MAX;
	for (int j = 0; j < samplePoints.size(); j++) {
		if (M[curOffset + j] < min) {
			sampleIndices[anchorPoints.size() - 1] = j;
		}
	}

	for (int i = anchorPoints.size() - 2; i >= 0; i--) {
		sampleIndices[i] = record[samplePoints.size()*(i + 1) + sampleIndices[i + 1]];
	}

	free(M);
	free(record);

	return sampleIndices;
}

void StructurePropagation::SetParm(int _blocksize, int _samplestep, int _iscurve) {
	this->blockSize = _blocksize;
	this->sampleStep = _samplestep;
	this->isCurve = _iscurve;
}

double StructurePropagation::calcE2(const Mat &mat, PointPos i1, PointPos i2) {
	int colLeft1, colLeft2, colRight1, colRight2;
	int rowUp1, rowUp2, rowDown1, rowDown2;
	Point p1 = pointManager.getPoint(i1);
	Point p2 = pointManager.getPoint(i2);
	if (p1.x > p2.x) {
		colLeft1 = 0;
		colLeft2 = p1.x - p2.x;
		colRight1 = blockSize - colLeft2;
		colRight2 = blockSize;
	}
	else {
		colLeft2 = 0;
		colLeft1 = p2.x - p1.x;
		colRight2 = blockSize - colLeft1;
		colRight1 = blockSize;
	}

	if (p1.y > p2.y) {
		rowUp1 = 0;
		rowUp2 = p1.y - p2.y;
		rowDown1 = blockSize - rowUp2;
		rowDown2 = blockSize;
	}
	else {
		rowUp2 = 0;
		rowUp1 = p2.y - p1.y;
		rowDown2 = blockSize - rowUp1;
		rowDown1 = blockSize;
	}
	if (colRight1 >= 0 && colRight2 >= 0 && rowDown1 >= 0 && rowDown2 >= 0) {
		double ssd = 0.0;
		int cols = colRight1 - colLeft1;
		int rows = rowDown1 - rowUp1;
		int xOffset1 = colLeft1 + p1.x - blockSize / 2;
		int xOffset2 = colLeft2 + p2.x - blockSize / 2;
		int yOffset1 = rowUp1 + p1.y - blockSize / 2;
		int yOffset2 = rowUp2 + p2.y - blockSize / 2;
		for (int i = 0; i < rows; i++) {
			const uchar *ptr1 = mat.ptr<uchar>(i + yOffset1);
			const uchar *ptr2 = mat.ptr<uchar>(i + yOffset2);
			for (int j = 0; j < cols; j++) {
				double diff = ptr1[j + xOffset1] - ptr2[j + xOffset2];
				ssd += diff*diff;
			}
		}
		return ssd / (cols * rows);
	}
	else {
		return 0;
	}
}

double StructurePropagation::calcEs(PointPos i, PointPos xi) {
	if (isCurve) {
		vector<Point> points1, points2;
		vector<double> minDistance1, minDistance2;
		pointManager.getPointsinPatch(i, points1);
		pointManager.getPointsinPatch(xi, points2);
		minDistance1.resize(points1.size());
		minDistance2.resize(points2.size());
		for (int i = 0; i < points1.size(); i++) {
			for (int j = 0; j < points2.size(); j++) {
				double diffx = points1[i].x - points2[j].x;
				double diffy = points1[i].y - points2[j].y;
				double distance = diffx*diffx + diffy*diffy;
				if (distance < minDistance1[i]) {
					minDistance1[i] = distance;
				}
				if (distance < minDistance2[j]) {
					minDistance2[j] = distance;
				}
			}
		}
		double es1 = 0.0, es2 = 0.0;
		for (int i = 0; i < minDistance1.size(); i++) {
			es1 += minDistance1[i];
		}
		for (int i = 0; i < minDistance2.size(); i++) {
			es2 += minDistance2[i];
		}
		return es1 / minDistance1.size() + es2 / minDistance2.size();
	}
	else {
		return 0.0;
	}
}

double StructurePropagation::calcEi(const Mat &mat, PointPos i, PointPos xi) {
	if (pointManager.nearBoundary(i)) {
		int offset1 = blockSize / 2;
		int offset2 = blockSize - offset1;
		double ssd = 0.0;
		Point pi = pointManager.getPoint(i);
		Point pxi = pointManager.getPoint(xi);
		for (int i = -offset1; i < offset2; i++) {
			const uchar *ptri = mat.ptr<uchar>(i + pi.y);
			const uchar *ptrxi = mat.ptr<uchar>(i + pxi.y);
			for (int j = -offset1; j < offset2; j++) {
				double diff = ptri[j + pi.x] - ptrxi[j + pxi.x];
				ssd += diff*diff;
			}
		}
		return ssd / (blockSize * blockSize);
	}
	else {
		return 0.0;
	}
}

