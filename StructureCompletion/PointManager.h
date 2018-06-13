#include "stdafx.h"
#include <hash_map>
#include <memory>

class Endpoints {
public:
	int trueLineIndex;
	int startIndex;
	int endIndex;
};

class PointPos {
public:
	const int lineIndex;
	const int pointIndex;
	PointPos(int lineIndex = -1, int pointIndex = -1): lineIndex(lineIndex), pointIndex(pointIndex) {
	}
};

class MElement {
public:
	double value;
	int xi;
	MElement() {

	}
	MElement(double value, int xi) {
		this->value = value;
		this->xi = xi;
	}
};

class Edge;
class Node {
private:
	list<Edge> edges;
	int edgeNum;
public:
	const PointPos p;
	const ushort id;
	static ushort totalNum;
	Node(PointPos p) : p(p), id(totalNum++) {
		edgeNum = 0;
	}
	void insertEdge(Edge *e) {
		edges.push_front(*e);
		edgeNum++;
	}
	void eraseEdge(list<Edge>::iterator itor) {
		edges.erase(itor);
		edges.push_back(*itor);
		edgeNum--;
	}
	void getEdges(list<Edge> &edges) {
		edges = this->edges;
	}
	int getEdgeNum() {
		return edgeNum;
	}
};

class Edge {
public:
	ushort ni;
	ushort nj;
	shared_ptr<MElement> Mij;
	shared_ptr<MElement> Mji;

	Edge(ushort ni, ushort nj, int size): ni(ni), nj(nj) {
		Mij = shared_ptr<MElement>((MElement*)malloc(size * sizeof(MElement)), free);
		Mji = shared_ptr<MElement>((MElement*)malloc(size * sizeof(MElement)), free);
	}

	Edge(int size) {
		Mij = shared_ptr<MElement>((MElement*)malloc(size * sizeof(MElement)), free);
		Mji = shared_ptr<MElement>((MElement*)malloc(size * sizeof(MElement)), free);
	}
};

class PointManager {
public:
	PointManager() {
	}
	void reset(const vector<vector<Point>> &linePoints, const Mat1b &mask, int blockSize);
	Point getPoint(PointPos p);
	bool nearBoundary(PointPos p);
	void getPointsinPatch(PointPos p, vector<Point> &ret);
	Point *getLinePtr(int i, int *length);
	int getLineNum();
	// void getSamplePoints(vector<Point> &samples);
	// void getAnchorPoints(vector<Point> &anchors);
	// void getIntersection(vector<list<PointPos>> &intersections);
	void constructBPMap(map<int, list<PointPos>> &intersectingMap);


private:
	vector<vector<Point>> linePoints;
	Mat1b mask;
	int blockSize;
	vector<Endpoints> lineEnds;
	set<PointPos> boundaryPoints;
	map<int, list<PointPos>> intersectingMap;
	vector<list<Node>> nodeListBucket;
	void constructMap();
	vector<list<Node>::iterator> nodeList;

	bool nearBoundary(const Point &p);
	int calcHashValue(int x, int y);
	Node *getBPNext();
};