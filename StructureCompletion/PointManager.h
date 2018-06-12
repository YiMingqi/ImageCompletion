#include "stdafx.h"

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
	PointPos();
	PointPos(int lineIndex, int pointIndex): lineIndex(lineIndex), pointIndex(pointIndex) {
	}
};

class MElement {
public:
	double value;
	int xi;
	MElement();
	MElement(double value, int xi) {
		this->value = value;
		this->xi = xi;
	}
};

class Edge;
class Node {
public:
	const PointPos p;
	list<Edge> edges;
	Node();
	Node(PointPos p) : p(p) {
	}
};

class Edge {
public:
	const list<Node>::iterator ni;
	const list<Node>::iterator nj;
	vector<MElement> Mij;
	vector<MElement> Mji;

	Edge(list<Node>::iterator ni, list<Node>::iterator nj, int size): ni(ni), nj(nj) {
		Mij.resize(size);
		Mji.resize(size);
	}
};

class PointManager {
public:
	PointManager() {
	}
	void reset(vector<vector<Point>> linePoints, const Mat1b &mask, int blockSize);
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

	bool nearBoundary(const Point &p);
	int calcHashValue(int x, int y);
};