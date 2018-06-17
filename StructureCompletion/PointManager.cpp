#include "stdafx.h"
#include "PointManager.h"

#define getLineIndex(x) (x-1) >> 24
#define getPointIndex(x) (x-1) << 8 >> 8
#define visit(l, p) ((l << 24) | (p + 1))

ushort Node::totalNum = 0;

bool operator<(const PointPos &p1, const PointPos &p2) {
	return (p1.lineIndex == p2.lineIndex) ? p1.pointIndex < p2.pointIndex : p1.lineIndex < p2.lineIndex;
}

bool operator==(const Edge &e1, const Edge &e2) {
	return (e1.ni == e2.ni) & (e1.nj == e2.nj);
}


void PointManager::reset(const vector<vector<Point>> &linePoints, const Mat1b &mask, int blockSize) {
	this->linePoints = linePoints;
	this->mask = mask;
	this->blockSize = blockSize;
	lineEnds.clear();
	boundaryPoints.clear();
	intersectingMap.clear();
	nodes.clear();
	propagationStack.clear();

	Mat visitMat = Mat::zeros(mask.rows, mask.cols, CV_32SC1);

	// record lines crossing the mask area and the patch center points near boundary
	bool inMask = false;
	Endpoints endpoints;
	for (int j = 0; j < linePoints.size(); j++) {
		int i;
		for (i = 0; i < linePoints[j].size(); i++) {
			int y = linePoints[j][i].y;
			int x = linePoints[j][i].x;
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
				if (nearBoundary(linePoints[j][i], false)) {
					boundaryPoints.insert(PointPos(lineEnds.size(), i));
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

bool PointManager::nearBoundary(const Point &p, bool isSample) {
	int leftBound = MAX(p.x - blockSize / 2, 0);
	int rightBound = MIN(p.x + blockSize - blockSize / 2, mask.cols);
	int upBound = MAX(p.y - blockSize / 2, 0);
	int downBound = MIN(p.y + blockSize - blockSize / 2, mask.rows);
	const uchar *upPtr = mask.ptr<uchar>(upBound);
	const uchar *downPtr = mask.ptr<uchar>(downBound - 1);
	for (int i = leftBound; i < rightBound; i++) {
		if (!upPtr[i] == isSample || !downPtr[i] == isSample) {
			return true;
		}
	}
	for (int i = upBound + 1; i < downBound - 1; i++) {
		if (!mask.at<uchar>(i, leftBound) == isSample || !mask.at<uchar>(i, rightBound - 1) == isSample) {
			return true;
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
	int hashValue = calcHashValue(center.x, center.y);
	list<PointPos> pointPositions;
	if (intersectingMap.count(hashValue)) {
		pointPositions = intersectingMap[hashValue];
	}
	else {
		pointPositions.push_back(p);
	}
	for (list<PointPos>::iterator p = pointPositions.begin(); p != pointPositions.end(); p++) {
		Endpoints endPoints = lineEnds[p->lineIndex];
		Point *points = &linePoints[endPoints.trueLineIndex][0];
		int beginIndex;
		for (int i = p->pointIndex; i >= 0; i--) {
			if (points[i].x < leftBound || points[i].y < upBound || points[i].x >= rightBound || points[i].y >= downBound) {
				beginIndex = i + 1;
				break;
			}
		}
		for (int i = beginIndex; i < p->pointIndex; i++) {
			ret.push_back(points[i]);
		}
		for (int i = p->pointIndex; i < linePoints[endPoints.trueLineIndex].size(); i++) {
			if (points[i].x < leftBound || points[i].y < upBound || points[i].x >= rightBound || points[i].y >= downBound) {
				break;
			}
			else {
				ret.push_back(points[i]);
			}
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

void PointManager::constructBPMap() {
	map<int, list<PointPos>>::iterator mapItor;
	list<shared_ptr<Node>> BFSstack;
	vector<vector<ushort>> pointVisitedMarks(linePoints.size());
	vector<list<shared_ptr<Node>>> nodeListBucket(4);

	// initialize the visit map
	pointVisitedMarks.resize(linePoints.size());
	for (int i = 0; i < linePoints.size(); i++) {
		pointVisitedMarks[i].resize(linePoints[i].size());
	}

	int total = 0;
	for (int i = 0; i < linePoints.size(); i++) {
		total += linePoints[i].size();
	}
	// reserve enough space for the node table
	nodes.reserve(total / blockSize);
	// skip the entry numbered with 0
	nodes.resize(1);

	// enqueue all the intersections
	vector<int> lineVisitedMarks(lineEnds.size());
	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++) {
		// enqueue the intersection (choose one point's position to represent all)
		BFSstack.push_back(make_shared<Node>(*(mapItor->second.begin())));
		// mark all intersecting points as visited
		list<PointPos>::iterator listItor = mapItor->second.begin();
		for (; listItor != mapItor->second.end(); listItor++) {
			pointVisitedMarks[lineEnds[listItor->lineIndex].trueLineIndex][listItor->pointIndex] = Node::totalNum;
			lineVisitedMarks[listItor->lineIndex] = 1;
		}
	}

	// enqueue a node in curves not intersecting with others
	for (int i = 0; i < lineEnds.size(); i++) {
		if (lineVisitedMarks[i] == 0) {
			// choose the first anchor point as starting point of propagation
			Endpoints endpoints = lineEnds[i];
			// enqueue the node
			BFSstack.push_back(make_shared<Node>(PointPos(i, endpoints.startIndex)));
			// mark the point as visited
			pointVisitedMarks[endpoints.trueLineIndex][endpoints.startIndex] = Node::totalNum;
		}
	}

	// enqueue all the neighbor nodes of intersections
	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++) {
		shared_ptr<Node> n = *BFSstack.begin();
		list<PointPos>::iterator listItor = mapItor->second.begin();
		int neighborNum = 0;
		for (; listItor != mapItor->second.end(); listItor++) {
			neighborNum += addNeighbor(*n, *listItor, pointVisitedMarks, BFSstack);
		}
		// Enlarge nodeListBucket if necessary
		if (neighborNum > nodeListBucket.size()) {
			nodeListBucket.resize(mapItor->second.size() * 2);
		}
		nodeListBucket[neighborNum - 1].push_front(n);
		nodes.push_back(nodeListBucket[neighborNum - 1].begin());
		BFSstack.pop_front();
	}

	// start propagation
	while (BFSstack.size()) {
		shared_ptr<Node> n = *BFSstack.begin();
		int neighborNum = addNeighbor(*n, n->p, pointVisitedMarks, BFSstack);
		nodeListBucket[neighborNum - 1].push_front(n);
		nodes.push_back(nodeListBucket[neighborNum - 1].begin());
		BFSstack.pop_front();
	}

	//generate the sequence for message sending
	while (nodeListBucket[0].size() > 0) {
		shared_ptr<Node> n = *nodeListBucket[0].begin();
		assert(n->getEdgeNum() == 1);
		list<shared_ptr<Edge>>::iterator eItor = n->getEdgeBegin();
		nodes[n->id] = propagationStack.insert(propagationStack.end(), n);
		nodeListBucket[0].pop_front();
		// degrade the adjacent node
		int id = (*eItor)->getAnother(n->id);
		n = *nodes[id];
		int edgeNum = n->getEdgeNum();
		n->eraseEdge(*eItor);
		nodeListBucket[edgeNum - 1].erase(nodes[id]);
		if (edgeNum > 1) {
			nodes[id] = nodeListBucket[edgeNum - 2].insert(nodeListBucket[edgeNum - 2].end(), n);
		}
		else {
			nodes[id] = propagationStack.insert(propagationStack.end(), n);
		}
	}
	int i = 0;
	i++;
}

/* void PointManager::constructBPMap() {
	map<int, list<PointPos>>::iterator mapItor;
	list<Node> BFSstack;
	vector<vector<ushort>> pointVisitedMarks;
	pointVisitedMarks.resize(linePoints.size());
	for (int i = 0; i < linePoints.size(); i++) {
		pointVisitedMarks[i].resize(linePoints[i].size());
	}
	nodeListBucket.resize(4);
	int total = 0;
	for (int i = 0; i < linePoints.size(); i++) {
		total += linePoints[i].size();
	}
	nodeList.reserve(total / blockSize);
	nodeList.resize(1);

	int *lineVisitedMarks = (int*)malloc(lineEnds.size() * sizeof(int));
	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++) {
		BFSstack.push_back(Node(*(mapItor->second.begin())));
		//预着色,用于处理两个不同位置的交点patch重叠的情况
		list<PointPos>::iterator listItor = mapItor->second.begin();
		for (; listItor != mapItor->second.end(); listItor++) {
			pointVisitedMarks[lineEnds[listItor->lineIndex].trueLineIndex][listItor->pointIndex] = Node::totalNum;
			lineVisitedMarks[listItor->lineIndex] = 1;
		}
		if (2 * mapItor->second.size() > nodeListBucket.size()) {
			nodeListBucket.resize(2 * mapItor->second.size());
		}
	}

	for (mapItor = intersectingMap.begin(); mapItor != intersectingMap.end(); mapItor++) {
		list<PointPos>::iterator listItor = mapItor->second.begin();
		for (; listItor != mapItor->second.end(); listItor++) {
			addNeighbor(*BFSstack.begin(), *listItor, pointVisitedMarks, BFSstack);
		}
		BFSstack.pop_front();
	}

	while (BFSstack.size()) {
		list<Node>::iterator itor = BFSstack.begin();
		addNeighbor(*itor, itor->p, pointVisitedMarks, BFSstack);
		BFSstack.pop_front();
	}

	Edge edge;
	for (int i = 0; i < lineEnds.size(); i++) {
		if (lineVisitedMarks[i] == 0) {
			Endpoints endPoints = lineEnds[i];
			int endIndex = endPoints.endIndex - blockSize / 2;
			for (int j = endPoints.startIndex; j < endIndex; j += blockSize / 2) {
				Node tmpNode = Node(PointPos(i, j));
				if (edge.ni > 0) {
					tmpNode.insertEdge(edge);
				}
				edge.ni = tmpNode.id;
				edge.nj = tmpNode.id + 1;
				tmpNode.insertEdge(edge);
				nodeListBucket[tmpNode.getEdgeNum() - 1].push_front(tmpNode);
			}
		}
	}

	free(lineVisitedMarks);

}*/

/*void PointManager::addNeighbor(Node &n, const PointPos &pos, vector<vector<ushort>> &visitedMark, list<Node> &BFSstack) {
	Endpoints endpoints = lineEnds[pos.lineIndex];
	int lineIndex = endpoints.trueLineIndex;
	int pointIndex = pos.pointIndex;
	int prePointIndex = pointIndex - blockSize / 2;
	int nextPointIndex = pointIndex + blockSize / 2;
	int bucketEntry = n.getEdgeNum() - 1;
	if (prePointIndex >= endpoints.startIndex) {
		int i;
		for (i = pointIndex - 1; i >= prePointIndex; i--) {
			if (visitedMark[lineIndex][i] && nodeList.size() > visitedMark[lineIndex][i]) {
				Edge tmpEdge(n.id, visitedMark[lineIndex][i]);
				n.insertEdge(tmpEdge);
				nodeList[visitedMark[lineIndex][i]]->insertEdge(tmpEdge);
				break;
			}
		}
		if (i == prePointIndex - 1) {
			BFSstack.push_back(Node(PointPos(lineIndex, prePointIndex)));
			visitedMark[lineIndex][i] = Node::totalNum;
		}
		bucketEntry++;
	}
	if (nextPointIndex < endpoints.endIndex) {
		int i;
		for (i = pointIndex + 1; i <= nextPointIndex; i++) {
			if (visitedMark[lineIndex][i] && nodeList.size() > visitedMark[lineIndex][i]) {
				Edge tmpEdge(n.id, visitedMark[lineIndex][i]);
				n.insertEdge(tmpEdge);
				nodeList[visitedMark[lineIndex][i]]->insertEdge(tmpEdge);
				break;
			}
		}
		if (i == nextPointIndex + 1) {
			BFSstack.push_back(Node(PointPos(lineIndex, nextPointIndex)));
			visitedMark[lineIndex][i] = Node::totalNum;
		}
		bucketEntry++;
	}
	nodeListBucket[bucketEntry].push_front(n);
	nodeList.push_back(nodeListBucket[bucketEntry].begin());
}*/


int PointManager::addNeighbor(Node &n, const PointPos &pos, vector<vector<ushort>> &visitedMark, list<shared_ptr<Node>> &BFSstack) {
	Endpoints endpoints = lineEnds[pos.lineIndex];
	int lineIndex = endpoints.trueLineIndex;
	int pointIndex = pos.pointIndex;
	int prePointIndex = pointIndex - blockSize / 2;
	int nextPointIndex = pointIndex + blockSize / 2;
	int neighborNum = 0;

	// check the point before current anchor point
	if (prePointIndex >= endpoints.startIndex) {
		int i;
		// try choosing a existed anchor point as its neighbor
		for (i = prePointIndex; i < pointIndex; i++) {
			if (visitedMark[lineIndex][i]) {
				if (nodes.size() > visitedMark[lineIndex][i]) {
					// add an edge between two points
					shared_ptr<Edge> tmpEdge = make_shared<Edge>(n.id, visitedMark[lineIndex][i]);
					n.push_front(tmpEdge);
					(*nodes[visitedMark[lineIndex][i]])->push_back(tmpEdge);
				}
				break;
			}
		}
		// no existed point can be chosen, construct a new anchor point and enqueue it
		if (i == pointIndex) {
			BFSstack.push_back(make_shared<Node>(PointPos(lineIndex, prePointIndex)));
			visitedMark[lineIndex][prePointIndex] = Node::totalNum;
		}
		neighborNum++;
	}

	// check the point behind current anchor point
	if (nextPointIndex < endpoints.endIndex) {
		int i;
		// try choosing a existed anchor point as its neighbor
		for (i = nextPointIndex; i > pointIndex; i--) {
			if (visitedMark[lineIndex][i]) {
				if (nodes.size() > visitedMark[lineIndex][i]) {
					// add an edge between two points
					shared_ptr<Edge> tmpEdge = make_shared<Edge>(n.id, visitedMark[lineIndex][i]);
					n.push_front(tmpEdge);
					(*nodes[visitedMark[lineIndex][i]])->push_back(tmpEdge);
				}
				break;
			}
		}
		// no existed point can be chosen, construct a new anchor point and enqueue it
		if (i == pointIndex) {
			BFSstack.push_back(make_shared<Node>(PointPos(lineIndex, nextPointIndex)));
			visitedMark[lineIndex][nextPointIndex] = Node::totalNum;
		}
		neighborNum++;
	}
	return neighborNum;
}

void PointManager::getPropstackItor(list<shared_ptr<Node>>::iterator &begin, list<shared_ptr<Node>>::iterator &end) {
	begin = propagationStack.begin();
	end = propagationStack.end();
}

void PointManager::getPropstackReverseItor(list<shared_ptr<Node>>::reverse_iterator &begin, list<shared_ptr<Node>>::reverse_iterator &end) {
	begin = list<shared_ptr<Node>>::reverse_iterator(propagationStack.end())++;
	end = list<shared_ptr<Node>>::reverse_iterator(propagationStack.begin());
}

unique_ptr<Node> PointManager::getBPNext() {
	/*if (nodeListBucket[0].size() > 0) {
		list<Node>::iterator itor = nodeListBucket[0].begin();
		assert(itor->getEdgeNum() == 1);
		list<Edge> e;
		itor->getEdges(e);
		list<Edge>::iterator eItor = e.begin();
		int id = (eItor->ni == itor->id) ? eItor->nj : eItor->ni;
		int edgeNum = nodeList[id]->getEdgeNum();
		nodeList[id]->eraseEdge(eItor);
		Node tmpNode = *nodeList[id];
		nodeListBucket[edgeNum - 1].erase(nodeList[id]);
		nodeListBucket[edgeNum - 2].push_front(tmpNode);
		nodeList[id] = nodeListBucket[edgeNum - 2].begin();
		nodeListBucket[0].pop_front();
		unique_ptr<Node> ret(new Node(*itor));
		return ret;
	}
	else {
		return NULL;
	}*/
	return NULL;
}

void PointManager::getSamplePoints(vector<PointPos> &samples, int sampleStep) {
	samples.clear();
	int lineIndex = 0;
	Endpoints endpoints = lineEnds[0];
	int total = 0;
	for (int i = 0; i < linePoints.size(); i++) {
		total += linePoints[i].size();
	}
	for (int i = 0; i < lineEnds.size(); i++) {
		total -= (lineEnds[i].endIndex - lineEnds[i].startIndex);
	}
	samples.reserve(total / sampleStep);
	for (int i = 0; i < linePoints.size(); i++) {
		int beginIndex = blockSize; //ensure all samples have complete line segments
		int endIndex;
		while (endpoints.trueLineIndex == i) {
			endIndex = endpoints.startIndex;
			for (int j = endIndex - 1; j >= beginIndex; j -= sampleStep) {
				if (j == beginIndex) {
					int c = 0;
					c++;
				}
				if (!nearBoundary(linePoints[i][j], true)) {
					samples.push_back(PointPos(lineIndex, j));
				}
			}
			beginIndex = endpoints.endIndex;
			++lineIndex;
			if (lineIndex >= lineEnds.size()) {
				break;
			}
			endpoints = lineEnds[lineIndex];
		}
		endIndex = linePoints[i].size() - blockSize; //ensure all samples have complete line segments
		for (int j = endIndex - 1; j >= beginIndex; j -= sampleStep) {
			if (!nearBoundary(linePoints[i][j], true)) {
				samples.push_back(PointPos(lineIndex - 1, j));
			}
		}
	}
	samples.shrink_to_fit();
}

void PointManager::getAnchorPoints(vector<PointPos> &anchors) {
	anchors.clear();
	int lineIndex = 0;
	Endpoints endpoints = lineEnds[0];
	int total = 0;
	for (int i = 0; i < lineEnds.size(); i++) {
		total += (lineEnds[i].endIndex - lineEnds[i].startIndex);
	}
	anchors.reserve(total / (blockSize / 2));
	for (int i = 0; i < linePoints.size(); i++) {
		while (endpoints.trueLineIndex == i) {
			for (int j = endpoints.startIndex; j < endpoints.endIndex; j += blockSize / 2) {
				anchors.push_back(PointPos(lineIndex, j));
			}
			++lineIndex;
			if (lineIndex >= lineEnds.size()) {
				break;
			}
			endpoints = lineEnds[lineIndex];
		}
	}
	anchors.shrink_to_fit();
}