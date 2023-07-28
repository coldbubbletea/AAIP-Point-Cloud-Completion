#pragma once
#include "OptAssign.h"

class COptUpdate :
    public COptAssign
{
public:
    COptUpdate(void);
    ~COptUpdate(void);

    void updateSIA();

    // update process
    void reformSubGraph(set<int>* previous_assign);
    void handleaffectedB(int &curdijkloop);
    int compDijkstra_rE(int &curloop, vector<int> &visited);
    int augmentFlow_rE(int lastid);
    void rebuildPsi(int &curdijkloop);
    void fixNewPsi(int &curdijkloop);
    void updatePotential_rE(int mindist, vector<int> &visited);

    // old method
    void updateSIA_slow();
    void reformSubGraph_slow(set<int>* previous_assign);
};
