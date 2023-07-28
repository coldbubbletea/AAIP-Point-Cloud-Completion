#pragma once
#include "OptUpdate.h"

//#define DIM 2

/**********************************/
/* Platform settings
/**********************************/
inline double gettime() {
    return (double)(clock() / (double)CLOCKS_PER_SEC);
}

struct Point2D { float v[3]; };
struct Rect2D { float v[6]; };

#define TIMESFACTOR 5931.0f
#define FloatToInt(a) (int)((double)(sqrt(a))*TIMESFACTOR)
#define IntToFloat(a) a/TIMESFACTOR

#define MAX_FLOAT numeric_limits<float>::max()
#define MAX_INT numeric_limits<int>::max()
#define MAX_UNSIGNEDINT numeric_limits<unsigned int>::max()

struct HeapEntry {
	int level,son;
	float key;
	//float bounces[DIM*2];
};

struct HeapComp {
	bool operator () (HeapEntry left,HeapEntry right) const
	{ return left.key > right.key; }
};

typedef	priority_queue<HeapEntry,deque<HeapEntry>,HeapComp> Heap;
typedef vector<HeapEntry> EntryList;

class RTree;

class CSpatialCCA : public COptUpdate
{
public:
	CSpatialCCA(void);
	~CSpatialCCA(void);

	// virtual functions
	void initLists();
	bool getNextList(int id, fHeap<int>::elem &el);

	// tool functions
	void readQry( int qrysize,int objsize,float* v);
	void testRTree();
	void readObj(RTree* RT,char* objfile, int objsize);
	void readObjbounces(char* objfile);
	void printstat(const char *processtype, char* objfile, int output=0);
	void readObjUpdates(char* objupdatefile);
    void readAssignment(char* assignfile);

	float MINDIST_nosqrt(Point2D &point, float *bounces, int dim=3);

private:
	vector<Point2D> Qry;
	vector<Point2D> Obj;
	vector<Heap> QryH;
	vector<vector<pair<int,double>>> bf_v;
    map<int,int> delid;

    int no_objupdates;
    int no_qryupdates;

	RTree *rt;
    int output;
};
