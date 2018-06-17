#include "stdafx.h"
#include "StructurePropagation.h"


void StructurePropagation::Run(const Mat1b &_mask, const Mat& _img, vector<vector<Point>> &linePoints, Mat& result) {

	// Calculate gray map
	Mat grayMat = Mat::zeros(_img.rows, _img.cols, CV_8UC1);
	for (int i = 0; i < _img.rows; i++) {
		for (int j = 0; j < _img.cols; j++) {
			Vec3b tmp = _img.at<Vec3b>(i, j);
			grayMat.at<uchar>(i, j) = (uchar)((114 * tmp[0] + 587 * tmp[1] + 299 * tmp[2] + 500) / 1000);
		}
	}

	// pointManager = PointManager();
	pointManager.reset(linePoints, grayMat, blockSize);
	vector<PointPos> samplePoints;
	pointManager.getSamplePoints(samplePoints, sampleStep);
	// No enough sample points or anchor points
	if (samplePoints.size() == 0){
		return;
	}

	int *sampleIndices;
	vector<PointPos> anchorPoints;
	sampleIndices = DP(samplePoints, anchorPoints, grayMat);
	// sampleIndices = BP(samplePoints, anchorPoints, grayMat);
	for (int i = 0; i < anchorPoints.size(); i++) {
		cout << sampleIndices[i] << ", ";
	}
	getResult(sampleIndices, samplePoints, anchorPoints, result);
}

void StructurePropagation::getResult(int *sampleIndices, const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, Mat& result) {
	// copy all sample patches to corresponding anchor pathces
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

int *StructurePropagation::BP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat) {
	// do initialization
	pointManager.constructBPMap();
	int size = pointManager.getPropstackSize();
	anchorPoints.reserve(size);

	list<shared_ptr<Node>>::iterator itor;
	list<shared_ptr<Node>>::iterator end;
	pointManager.getPropstackItor(itor, end);
	// receive message sent from other neighbors
	for (; itor != end; itor++) {
		shared_ptr<Node> n = *itor;
		// calculate message for next neighbor (the node that enqueued this node)
		calcMij(*n, n->getEdgeBegin(), mat, samplePoints);
	}

	// send updated message back to neighbors
	int *sampleIndices = (int*)malloc(size * sizeof(int));
	double *cur = (double*)malloc(samplePoints.size() * sizeof(double));

	list<shared_ptr<Node>>::reverse_iterator rev_itor;
	list<shared_ptr<Node>>::reverse_iterator rev_end;
	pointManager.getPropstackReverseItor(rev_itor, rev_end);
	for (int i = 0; rev_itor != rev_end; rev_itor++, i++) {
		shared_ptr<Node> n = *rev_itor;
		list<shared_ptr<Edge>>::iterator begin = n->getEdgeBegin();
		list<shared_ptr<Edge>>::iterator end = n->getEdgeEnd();
		list<shared_ptr<Edge>>::iterator itor = begin;

		anchorPoints.push_back(n->p);
		// calculate all the send-back message
		for (itor++; itor != end; itor++) {
			calcMij(*n, itor, mat, samplePoints);
		}
		int minIndex;
		double min = INT_MAX;
		// calculate E1 for all possible xi
		for (int i = 0; i < samplePoints.size(); i++) {
			cur[i] = ks * calcEs(n->p, samplePoints[i]) + ki * calcEi(mat, n->p, samplePoints[i]);
		}
		// add up all messages sent to this node
		for (itor = begin; itor == end; itor++) {
			double **toMptr = (*itor)->getMbyTo(n->id);
			for (int i = 0; i < samplePoints.size(); i++) {
				cur[i] += (*toMptr)[i];
			}
		}
		// find out the optimal xi
		for (int i = 0; i < samplePoints.size(); i++) {
			if (cur[i] < min) {
				min = cur[i];
				minIndex = i;
			}
		}
		sampleIndices[i] = minIndex;
	}

	// release resources
	pointManager.getPropstackItor(itor, end);
	for (; itor != end; itor++) {
		shared_ptr<Node> n = *itor;
		list<shared_ptr<Edge>>::iterator edgeItor = n->getEdgeBegin();
		list<shared_ptr<Edge>>::iterator end = n->getEdgeEnd();
		for (; edgeItor != end; edgeItor++) {
			double **M = (*edgeItor)->getMbyFrom(n->id);
			if (*M != NULL) {
				free(*M);
			}
		}
	}
	free(cur);

	return sampleIndices;
}

void StructurePropagation::calcMij(Node &n, const list<shared_ptr<Edge>>::iterator &edgeItor, const Mat &mat, const vector<PointPos> &samplePoints) {
	double **Mptr = (*edgeItor)->getMbyFrom(n.id);
	if (*Mptr == NULL) {
		*Mptr = (double*)malloc(samplePoints.size() * sizeof(double));
		for (int i = 0; i < samplePoints.size(); i++) {
			// calcalate Ei beforehand
			double E1 = ks * calcEs(n.p, samplePoints[i]) + ki * calcEi(mat, n.p, samplePoints[i]);
			PointPos tmpPos = pointManager.getPointPos((*edgeItor)->getAnother(n.id));
			for (int j = 0; j < samplePoints.size(); j++) {
				// try updating tne minimal value of each item in Mij
				double E2 = calcE2(mat, tmpPos, n.p, samplePoints[i], samplePoints[j]);
				if (E1 + E2 < (*Mptr)[j]) {
					(*Mptr)[j] = E1 + E2;
				}
			}
		}
		// add up the message sent from Mki (k != j)
		list<shared_ptr<Edge>>::iterator itor = n.getEdgeBegin();
		list<shared_ptr<Edge>>::iterator end = n.getEdgeEnd();
		for (; itor != end; itor++) {
			if (itor != edgeItor) {
				double **toMptr = (*itor)->getMbyTo(n.id);
				if (*toMptr == NULL) {
					assert(0);
				}
				for (int i = 0; i < samplePoints.size(); i++) {
					(*Mptr)[i] += (*toMptr)[i];
				}
			}
		}
	}
}

int *StructurePropagation::DP(const vector<PointPos> &samplePoints, vector<PointPos> &anchorPoints, const Mat &mat) {
	pointManager.getAnchorPoints(anchorPoints);

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
		// calculate the min value for each xi
		// xi = j
		for (int j = 0; j < samplePoints.size(); j++) {
			double E1 = ks * calcEs(anchorPoints[i], samplePoints[j]) +
				ki * calcEi(mat, anchorPoints[i], samplePoints[j]);
			double min = INT_MAX;
			// choose optimal x(i-1)
			// x(i-1) = k
			for (int k = 0; k < samplePoints.size(); k++) {
				double tmp = calcE2(mat, anchorPoints[i], anchorPoints[i - 1], samplePoints[j], samplePoints[k]) + M[preOffset + k];
				if (tmp < min) {
					record[samplePoints.size()*i + j] = k;
					min = tmp;
				}
			}
			M[curOffset + j] = E1 + min;
		}
	}

	// find out the optimal xi of last anchor point
	int *sampleIndices = (int*)malloc(anchorPoints.size() * sizeof(int));
	double min = INT_MAX;
	for (int j = 0; j < samplePoints.size(); j++) {
		if (M[curOffset + j] < min) {
			sampleIndices[anchorPoints.size() - 1] = j;
		}
	}

	// trace back
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

double StructurePropagation::calcE2(const Mat &mat, const PointPos &i1, const PointPos &i2, const PointPos &xi1, const PointPos &xi2) {
	int colLeft1, colLeft2, colRight1, colRight2;
	int rowUp1, rowUp2, rowDown1, rowDown2;
	Point p1 = pointManager.getPoint(i1);
	Point p2 = pointManager.getPoint(i2);
	Point px1 = pointManager.getPoint(xi1);
	Point px2 = pointManager.getPoint(xi2);
	// calculate the relative offsets of overlapping area's bounderies
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
		// calculate the absolute cooordinates of boundaries
		int xOffset1 = colLeft1 + px1.x - blockSize / 2;
		int xOffset2 = colLeft2 + px2.x - blockSize / 2;
		int yOffset1 = rowUp1 + px1.y - blockSize / 2;
		int yOffset2 = rowUp2 + px2.y - blockSize / 2;
		for (int i = 0; i < rows; i++) {
			const uchar *ptr1 = mat.ptr<uchar>(i + yOffset1);
			const uchar *ptr2 = mat.ptr<uchar>(i + yOffset2);
			for (int j = 0; j < cols; j++) {
				double diff = ptr1[j + xOffset1] - ptr2[j + xOffset2];
				ssd += diff*diff;
			}
		}
		// do normlization
		if (ssd != 0) {
			ssd += 0.0;
		}
		return ssd / (cols * rows);
	}
	else {
		// no overlapping part
		return 0.0;
	}
}

double StructurePropagation::calcEs(const PointPos &i, const PointPos &xi) {
	if (isCurve) {
		vector<Point> points1(30), points2(30);
		vector<int> minDistance1, minDistance2;
		Point pi = pointManager.getPoint(i);
		Point pxi = pointManager.getPoint(xi);
		int offsetx = pxi.x - pi.x;
		int offsety = pxi.y - pi.y;
		// get points of curve segment contained in patch
		pointManager.getPointsinPatch(i, points1);
		pointManager.getPointsinPatch(xi, points2);
		minDistance1.resize(points1.size());
		minDistance2.resize(points2.size());
		// initialize minimal distance
		for (int i = 0; i < points1.size(); i++) {
			minDistance1[i] = INT_MAX;
		}
		for (int i = 0; i < points2.size(); i++) {
			minDistance2[i] = INT_MAX;
		}
		// calculate the minimal distances for points in curve segments
		for (int i = 0; i < points1.size(); i++) {
			for (int j = 0; j < points2.size(); j++) {
				int diffx = points1[i].x - points2[j].x + offsetx;
				int diffy = points1[i].y - points2[j].y + offsety;
				int distance = diffx*diffx + diffy*diffy;
				if (distance < minDistance1[i]) {
					minDistance1[i] = distance;
				}
				if (distance < minDistance2[j]) {
					minDistance2[j] = distance;
				}
			}
		}
		int es1 = 0, es2 = 0;
		for (int i = 0; i < minDistance1.size(); i++) {
			es1 += minDistance1[i];
		}
		for (int i = 0; i < minDistance2.size(); i++) {
			es2 += minDistance2[i];
		}
		return (double)es1 / minDistance1.size() + (double)es2 / minDistance2.size();
	}
	else {
		// line segments in patches are always the same
		return 0.0;
	}
}

double StructurePropagation::calcEi(const Mat &mat, const PointPos &i, const PointPos &xi) {
	if (pointManager.nearBoundary(i)) {
		int offset1 = blockSize / 2;
		int offset2 = blockSize - offset1;
		int ssd = 0;
		int overlappingPixelNum = 0;
		Point pi = pointManager.getPoint(i);
		Point pxi = pointManager.getPoint(xi);
		for (int i = -offset1; i < offset2; i++) {
			const uchar *ptri = mat.ptr<uchar>(i + pi.y);
			const uchar *ptrxi = mat.ptr<uchar>(i + pxi.y);
			for (int j = -offset1; j < offset2; j++) {
				// filter out invalid pixels located in masked area
				if (ptri[j + pi.x] != 0) {
					int diff = ptri[j + pi.x] - ptrxi[j + pxi.x];
					overlappingPixelNum++;
					ssd += diff*diff;
				}
			}
		}
		// do nomalization
		return (double)ssd / overlappingPixelNum;
	}
	else {
		// target patch is contained in mask area totally
		return 0.0;
	}
}

