#pragma once

#include <iostream>
#include <fstream>
#include <string.h>

#include <vector>
#include <queue>
#include <list>
#include <limits>
#include <math.h>
#include <bitset>
#include <map>
#include <algorithm>
#include <bitset>
#include <set>

#include <time.h>
#include <math.h>

#include "nheap.h"

using namespace std;

#define DEBUG
#define TIMESFACTOR 5931.0f
#define FloatToInt(a) (int)((double)(sqrt(a)) * TIMESFACTOR)
#define IntToFloat(a) a / TIMESFACTOR
enum RunType
{
    SIA = 0,
    SIAUPDATE
};
enum DijkHeapUpdateType
{
    NEWNODE = 0,
    UPDATE,
    NOUPDATE
};
enum UpdateTags
{
    A,
    B,
    UP,
    DEL,
    INS
};

struct edge
{
    int fromid;
    int toid;
    int weight;
    int flow;
    int maxflow;
    int node_eidx;
    int node_reidx; // reverse
    // int maxflow;
};

struct node
{
    int psi;
    int flow;
    int maxflow;

    // for incremenal process
    fHeap<int>::elem best_el;

    // for dijkstra
    int mindist;
    int mineid;
    int augflow;
    int loop;

    // for re-dijkstra
    bool isUpdated;

    // for update process
    bool isUpdatedfromPrevious;

    // for update
    bool isCheck;

    vector<int> E;
    vector<int> rE; // we do not need rE in SIA

    map<int, int> preAssign;
};

class UpdateInfo
{
public:
    int update;
    int type;
    int id;
};

struct AssignE
{
    int from;
    int to;
    int cost;
    int cap;
    AssignE(int _from, int _to, int _cost, int _cap) : from(_from), to(_to), cost(_cost), cap(_cap) {}
};
class COptAssign
{

public:
    long long batch_size;
    long long item_size;
    vector<int> max_candidate_price;
    long long shutdown_time_milliseconds;
    COptAssign(void){};
    ~COptAssign(void);

    void initSIA(int runtype, int noA, int noB, int ACap[], int BCap[], int EdgeCap = 1);
    int runSIA(int * is_conflict=nullptr,int *curloopsent = nullptr, int curdijkloop = 0);
    void SIA_opt1(int fromid, vector<node>::iterator fromit);

    // tools
    void initDijkstra(int nid, int &curloop, vector<int> &visited);
    int compDijkstra_INC(int &curloop, fHeap<int> &globalH, vector<int> &visited);
    int updateDijkstra(int &curloop, fHeap<int> &globalH, vector<int> &visited);
    int updateMinDist(int curloop, int eid, vector<edge>::iterator eit,
                      vector<node>::iterator fromit, vector<node>::iterator toit,
                      int augflow, int cost);
    int insertEdgefromHeap(int curloop, fHeap<int> &globalH);
    int augmentFlow(int lastid, int &cost);
    void updatePotential(int mindist, vector<int> &visited);
    void updateHeaps(int curloop, int toid, vector<node>::iterator toit, int ret);
    void insertEdge(int fromid, int toid, int dist, vector<node>::iterator fromit, vector<node>::iterator toit);

    void insert_eid(vector<int> &eids, int &node_eidx, int eid);
    void remove_eid(vector<int> &eids, int del_eidx, int &updated_eidx);

    void insertAssign(int fromid, int toid, int dist, vector<node>::iterator fromit, vector<node>::iterator toit);
    void removeAssign(int eid, int toid);
    long long IPA(int *assignment, int *assignment_inv, float *price, int qrysize);

    int removeFeasibleAssign(int a, int b, int delta);
    // constraint assignment
    void computeConstraint();

    // inline functions
    inline int get_toid(int nid, vector<edge>::iterator eit)
    {
        return (eit->fromid == nid) ? eit->toid : eit->fromid;
    }
    inline int getNoAugment(vector<edge>::iterator eit,
                            vector<node>::iterator fromit, vector<node>::iterator toit)
    {
        // int flow=min(edgecap, min(fromit->maxflow,-toit->maxflow))-eit->flow; // THE bug!!!
        int flow = min(fromit->maxflow, -toit->maxflow) - eit->flow;
        if (fromit->maxflow > 0 && flow > 0)
            return flow;
        else if (fromit->maxflow < 0 && eit->flow > 0)
            return eit->flow;
        else
            return -1;
    }
    inline int calCost(vector<edge>::iterator eit, vector<node>::iterator fromit, vector<node>::iterator toit)
    {
        return (fromit->maxflow > 0) ? eit->weight - fromit->psi + toit->psi : -eit->weight - fromit->psi + toit->psi;
    }

    // virtual functions
    virtual void initLists(){};
    virtual bool getNextList(int id, fHeap<int>::elem &el) { return true; };
    virtual void initLists_r(){};
    virtual bool getNextList_r(int id, fHeap<int>::elem &el) { return true; };

protected:
    int noA, noB;
    // int Acap, Bcap, edgecap;
    int edgecap;
    int totalflows;
    int maxpsi;

    vector<node> nodes;
    vector<edge> edges;

    fHeap<int> dijkH;
    fHeap<int> updateH;

    vector<UpdateInfo> updates;

    int runtype;

    vector<AssignE> assignEdges;
    vector<vector<int>> asgn;
    vector<int> asgnHead;

    int currentcost;
    int remaincost;
    int currentflow;
    int feasiblecost;

    // double** costs;
    int **costs;
    int **sortedIdx;

    double *w1, *w2;
};
