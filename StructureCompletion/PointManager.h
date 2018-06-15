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

/*class MElement {
public:
	float value;
	int xi;
	MElement() {
		value = 0;
		xi = -1;
	}
	MElement(float value, int xi) {
		this->value = value;
		this->xi = xi;
	}
};*/

class Edge;
class Node {
private:
	list<Edge> edges;
	int edgeNum;
public:
	const PointPos p;
	const ushort id;
	static ushort totalNum;
	Node(PointPos p) : p(p), id(++totalNum) {
		edgeNum = 0;
	}
	void insertEdge(const Edge &e) {
		edges.push_front(e);
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
private:
	double *Mij;
	double *Mji;
public:
	ushort ni;
	ushort nj;

	Edge(ushort ni, ushort nj): ni(ni), nj(nj) {
	}

	Edge() {
	}

	inline ushort getAnother(ushort n) {
		if (n == ni) {
			return nj;
		}
		else {
			return ni;
		}
	}

	inline double **getMbyFrom(ushort from) {
		if (from == ni) {
			return &Mij;
		}
		else {
			return &Mji;
		}

	}

	inline double **getMbyTo(ushort to) {
		if (to == ni) {
			return &Mji;
		}
		else  {
			return &Mij;
		}
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
	void getSamplePoints(vector<PointPos> &samples, int sampleStep);
	void constructBPMap();
	unique_ptr<Node> getBPNext();
	void getAnchorPoints(vector<PointPos> &anchors);


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
	void addNeighbor(Node &n, const PointPos &pos, vector<vector<ushort>> &visitedMark, list<Node> &BFSstack);
};