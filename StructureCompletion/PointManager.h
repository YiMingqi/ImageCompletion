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

class Edge;
class Node {
public:
	const PointPos p;
	const ushort id;
	list<shared_ptr<Edge>> edges;
	static ushort totalNum;
	Node(PointPos p) : p(p), id(++totalNum) {
	}
	void getEdges(list<shared_ptr<Edge>> &edges) {
		edges = this->edges;
	}
};

class Edge {
private:
	double *Mij;
	double *Mji;
public:
	shared_ptr<Node> ni;
	shared_ptr<Node> nj;

	Edge(shared_ptr<Node> ni, shared_ptr<Node> nj) : ni(ni), nj(nj) {
		Mij = Mji = NULL;
	}

	Edge() {
		Mij = Mji = NULL;
	}

	inline shared_ptr<Node> getAnother(shared_ptr<Node> n) {
		if (n == ni) {
			return nj;
		}
		else {
			return ni;
		}
	}

	inline double **getMbyFrom(shared_ptr<Node> from) {
		if (from == ni) {
			return &Mij;
		}
		else {
			return &Mji;
		}

	}

	inline double **getMbyTo(shared_ptr<Node> to) {
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
	void reset(const vector<vector<Point>> &linePoints, const Mat1b &mask, int blockSize, set<shared_ptr<list<int>>> &lineSets);
	Point getPoint(PointPos p);
	bool nearBoundary(PointPos p);
	void getPointsinPatch(PointPos p, vector<Point> &ret);
	void getPointsinPatch(const PointPos &p, list<Point*> &begin, list<int> &length);
	void getSamplePoints(vector<PointPos> &samples, int sampleStep, list<int> &line);
	void getSamplePoints(vector<PointPos> &samples, int sampleStep);
	void constructBPMap(list<int> &line);
	void getAnchorPoints(vector<PointPos> &anchors, list<int> &line);
	void getNodesIterator(vector<shared_ptr<Node>>::iterator &begin, vector<shared_ptr<Node>>::iterator &end) {
		begin = ++nodes.begin();
		end = nodes.end();
	}
	shared_ptr<Node> getNode(ushort id) {
		return nodes[id];
	}

private:
	vector<vector<Point>> linePoints; //��¼�û����Ƶĵ����Ϣ
	Mat1b mask;
	int blockSize;
	vector<Endpoints> lineEnds; //���ڼ�¼����PointManager�ٴλ��ֺ���߶ε���β��Ϣ
	set<PointPos> boundaryPoints; //���ڼ�¼����patch��߽��ص���ê��
	map<int, list<PointPos>> intersectingMap; //���ڼ�¼���㣬��ֵΪ���ݽ������ʵ����������hashֵ
	vector<shared_ptr<Node>> nodes; //һ�Ÿ���Node id����node�ı���¼��Node�����ָ��
	// list<shared_ptr<Node>> propagationStack; //��¼BP�㷨����Ϣ���ݵ�˳��

	bool nearBoundary(const Point &p, bool isSample);
	int calcHashValue(int x, int y);
	void addNeighbor(shared_ptr<Node> &n, const PointPos &pos, vector<vector<ushort>> &visitedMark, list<Node> &BFSstack);
	int addNeighbor(shared_ptr<Node> &n, const PointPos &pos, vector<vector<ushort>> &visitedMark, list<shared_ptr<Node>> &neighbors);
};