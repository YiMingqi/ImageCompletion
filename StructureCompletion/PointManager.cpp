#include "stdafx.h"
#include "PointManager.h"

#define getLineIndex(x) (x-1) >> 24
#define getPointIndex(x) (x-1) << 8 >> 8
#define visit(l, p) ((l << 24) | (p + 1))

bool operator<(const PointPos &p1, const PointPos &p2) {
	return (p1.lineIndex == p2.lineIndex) ? p1.pointIndex < p2.pointIndex : p1.lineIndex < p2.lineIndex;
}

bool operator==(const Edge &e1, const Edge &e2) {
	return (e1.ni == e2.ni) & (e1.nj == e2.nj);
}

ushort Node::totalNum = 0;

void PointManager::reset(const vector<vector<Point>> &linePoints, const Mat1b &mask, int blockSize) {
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

struct {
	Node node;
	list<PointPos>::iterator itor;
};

void PointManager::constructBPMap(map<int, list<PointPos>> &intersectingMap) {
	map<int, list<PointPos>>::iterator mapItor;
	list<Node> BFSstack;
	vector<vector<ushort>> pointVisitedMarks;
	pointVisitedMarks.resize(linePoints.size());
	for (int i = 0; i < linePoints.size(); i++) {
		pointVisitedMarks[i].resize(linePoints[i].size());
		memset(&(pointVisitedMarks[i][0]), 0, linePoints[i].size() * sizeof(ushort));
	}
	nodeListBucket.resize(4);
	int total = 0;
	for (int i = 0; i < linePoints.size(); i++) {
		total += linePoints[i].size();
	}
	nodeList.reserve(total / blockSize);
	nodeList.resize(1);

	vector<Node> intersections;
	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++) {
		intersections.push_back(Node(*(mapItor->second.begin())));
		//预着色,用于处理两个不同位置的交点patch重叠的情况
		vector<ushort> marks = pointVisitedMarks[lineEnds[mapItor->second.begin()->lineIndex].trueLineIndex];
		int beginIndex = mapItor->second.begin()->pointIndex;
		for (int i = beginIndex - blockSize / 2 + 1; i < beginIndex + blockSize / 2; i++) {
			marks[i] = Node::totalNum;
		}
		if (2 * mapItor->second.size() > nodeListBucket.size()) {
			nodeListBucket.resize(2 * mapItor->second.size());
		}
	}
	int i = 0;
	//先用一个数组记录一下某条线段有没有被标记，没有被标记的线段从开头开始扩散，一直到结尾结束
	int *lineVisitedMarks = (int*)malloc(lineEnds.size() * sizeof(int));
	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++, i++) {
		list<PointPos>::iterator listItor = mapItor->second.begin();
		for (; listItor != mapItor->second.end(); listItor++) {
			addNeighbor(intersections[i], *listItor, pointVisitedMarks, BFSstack);
			lineVisitedMarks[listItor->lineIndex] = 1;
		}
	}

	while (BFSstack.size()) {
		list<Node>::iterator itor = BFSstack.begin();
		addNeighbor(*itor, itor->p, pointVisitedMarks, BFSstack);
		BFSstack.pop_front();
	}

	Edge edge(blockSize);
	for (int i = 0; i < lineEnds.size(); i++) {
		if (lineVisitedMarks[i] == 0) {
			Endpoints endPoints = lineEnds[i];
			vector<Point> points = linePoints[endPoints.trueLineIndex];
			int endIndex = endPoints.endIndex - blockSize / 2;
			for (int j = endPoints.startIndex; j < endIndex; j += blockSize / 2) {
				Node tmpNode = Node(PointPos(i, j));
				if (edge.ni >= 0) {
					tmpNode.insertEdge(edge);
				}
				edge.ni = tmpNode.id;
				edge.nj = tmpNode.id + 1;
				tmpNode.insertEdge(edge);
			}
		}
	}

	free(lineVisitedMarks);

}

void PointManager::addNeighbor(Node &n, const PointPos &pos, const vector<vector<ushort>> &visitedMark, list<Node> &BFSstack) {
	Endpoints endpoints = lineEnds[pos.lineIndex];
	int lineIndex = endpoints.trueLineIndex;
	int pointIndex = pos.pointIndex;
	int prePointIndex = pointIndex - blockSize / 2;
	int nextPointIndex = pointIndex + blockSize / 2;
	vector<ushort> marks = visitedMark[lineIndex];
	if (prePointIndex >= endpoints.startIndex) {
		ushort preNodeId = marks[prePointIndex];
		if (preNodeId && nodeList.size() > preNodeId) {
			Edge tmpEdge(n.id, preNodeId, blockSize);
			n.insertEdge(tmpEdge);
			nodeList[preNodeId]->insertEdge(tmpEdge);
		}
		else {
			BFSstack.push_back(Node(PointPos(lineIndex, prePointIndex)));
			for (int i = pointIndex - 1; i >= prePointIndex; i--) {
				marks[i] = Node::totalNum;
			}
		}
	}
	if (nextPointIndex < endpoints.endIndex) {
		ushort nextNodeId = marks[nextPointIndex];
		if (nextNodeId && nodeList.size() > nextNodeId) {
			Edge tmpEdge(n.id, nextNodeId, blockSize);
			n.insertEdge(tmpEdge);
			nodeList[nextNodeId]->insertEdge(tmpEdge);
		}
		else {
			BFSstack.push_back(Node(PointPos(lineIndex, nextPointIndex)));
			for (int i = pointIndex + 1; i <= nextPointIndex; i++) {
				marks[i] = Node::totalNum;
			}
		}
	}
	list<Node>::iterator itor;
	int bucketEntry = n.getEdgeNum() - 1;
	nodeListBucket[bucketEntry].push_front(n);
	nodeList.push_back(nodeListBucket[bucketEntry].begin());
}

Node *PointManager::getBPNext() {
	list<Node>::iterator itor = nodeListBucket[0].begin();
	assert(itor->getEdgeNum() == 1);
	list<Edge> e;
	itor->getEdges(e);
	list<Edge>::iterator eItor = e.begin();
	int id = (eItor->ni == itor->id) ? eItor->nj : eItor->ni;
	int edgeNum = nodeList[id]->getEdgeNum();
	nodeListBucket[edgeNum - 1].erase(nodeList[id]);
	nodeList[id]->eraseEdge(eItor);
	nodeListBucket[edgeNum - 2].push_front(*nodeList[id]);
	nodeList[id] = nodeListBucket[edgeNum - 2].begin();
	nodeListBucket[0].pop_front();
	return &(*itor);
}