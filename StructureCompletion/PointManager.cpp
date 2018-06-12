#include "stdafx.h"
#include "PointManager.h"

#define getLineIndex(x) (x-1) >> 24
#define getPointIndex(x) (x-1) << 8 >> 8
#define visit(l, p) ((l << 24) | (p+1))

bool operator<(const PointPos &p1, const PointPos &p2) {
	return (p1.lineIndex == p2.lineIndex) ? p1.pointIndex < p2.pointIndex : p1.lineIndex < p2.lineIndex;
}

void PointManager::reset(vector<vector<Point>> linePoints, const Mat1b &mask, int blockSize) {
	this->linePoints = linePoints;
	this->mask = mask;
	this->blockSize = blockSize;
	lineEnds.clear();
	boundaryPoints.clear();
	intersectingMap.clear();
	nodeListBucket.clear();

	Mat visitMat = Mat::zeros(mask.rows, mask.cols, CV_32SC1);

	// record lines crossing the mask area and the patch center points near boundary
	bool inMask = false;
	Endpoints endpoints;
	for (int j = 0; j < linePoints.size(); j++) {
		vector<Point> points = linePoints[j];
		int i;
		for (i = 0; i < points.size(); i++) {
			int y = points[i].y;
			int x = points[i].x;
			if (y < 0 || y >= mask.rows || x < 0 || x >= mask.cols) {
				continue;
			}
			else if (mask.at<uchar>(y, x)) {
				if (inMask == true) {
					endpoints.endIndex = i;
					lineEnds.push_back(endpoints);
					inMask = false;
				}
			}
			else {
				if (nearBoundary(points[i])) {
					boundaryPoints.insert(PointPos(j, i));
				}
				if (inMask == false) {
					endpoints.startIndex = i;
					endpoints.trueLineIndex = j;
					inMask = true;
				}
				int visitRecord = visitMat.at<int>(y, x);
				int lineIndex = getLineIndex(visitRecord);
				int pointIndex = getPointIndex(visitRecord);
				if (visitRecord != 0 && lineIndex != lineEnds.size()) {
					list<PointPos> *intersectingList = &intersectingMap[calcHashValue(x, y)];
					if (intersectingList->size() == 0) {
						intersectingList->push_back(PointPos(lineIndex, pointIndex));
					}
					intersectingList->push_back(PointPos(lineEnds.size(), i));
				}
				else {
					visitMat.at<int>(y, x) = visit(lineEnds.size(), i);
				}
			}
		}
		if (inMask == true) {
			inMask = false;
			endpoints.endIndex = i;
			lineEnds.push_back(endpoints);
		}
	}
}

bool PointManager::nearBoundary(const Point &p) {
	int leftBound = MAX(p.x - blockSize / 2, 0);
	int rightBound = MIN(p.x + blockSize - blockSize / 2, mask.cols);
	int upBound = MAX(p.y - blockSize / 2, 0);
	int downBound = MIN(p.y + blockSize - blockSize / 2, mask.rows);
	for (int i = upBound; i < downBound; i++) {
		const uchar *ptr = mask.ptr<uchar>(i);
		for (int j = leftBound; j < rightBound; j++) {
			if (ptr[j] != 0) {
				return true;
			}
		}
	}
	return false;
}

inline int PointManager::calcHashValue(int x, int y) {
	return x + y * mask.cols;
}

inline Point PointManager::getPoint(PointPos p) {
	return linePoints[lineEnds[p.lineIndex].trueLineIndex][p.pointIndex];
}

bool PointManager::nearBoundary(PointPos p) {
	return boundaryPoints.count(p);
}

void PointManager::getPointsinPatch(PointPos p, vector<Point> &ret) {
	Point center = getPoint(p);
	int leftBound = MAX(center.x - blockSize / 2, 0);
	int rightBound = MIN(center.x + blockSize - blockSize / 2, mask.cols);
	int upBound = MAX(center.y - blockSize / 2, 0);
	int downBound = MIN(center.y + blockSize - blockSize / 2, mask.rows);
	Endpoints endPoints = lineEnds[p.lineIndex];
	vector<Point> points = linePoints[endPoints.trueLineIndex];
	int beginIndex;
	for (int i = p.pointIndex; i >= endPoints.startIndex; i--) {
		Point point = points[i];
		if (point.x < leftBound || point.y < upBound || point.x >= rightBound || point.y >= downBound) {
			beginIndex = i + 1;
			break;
		}
	}
	for (int i = beginIndex; i < p.pointIndex; i++) {
		ret.push_back(points[i]);
	}
	for (int i = p.pointIndex; i < endPoints.endIndex; i++) {
		if (points[i].x < leftBound || points[i].y < upBound || points[i].x >= rightBound || points[i].y >= downBound) {
			break;
		}
		else {
			ret.push_back(points[i]);
		}
	}
}

int PointManager::getLineNum() {
	return lineEnds.size();
}

Point *PointManager::getLinePtr(int i, int *length) {
	if (i >= lineEnds.size()) {
		return NULL;
	}
	Endpoints endPoints = lineEnds[i];
	*length = endPoints.endIndex - endPoints.startIndex;
	return &linePoints[endPoints.trueLineIndex][endPoints.startIndex];
}

void constructBPMap(map<int, list<PointPos>> &intersectingMap) {
	map<int, list<PointPos>>::iterator itor;
	for (itor = intersectingMap.begin(); itor != intersectingMap.end(); itor++) {
		int lineIndex = getLineIndex((*itor).first);
		int pointIndex = getPointIndex((*itor).first);
		Node(PointPos(lineIndex, pointIndex));
	}
}