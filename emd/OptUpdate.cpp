#include "OptUpdate.h"

COptUpdate::COptUpdate(void)
{
}

COptUpdate::~COptUpdate(void)
{
}


void COptUpdate::reformSubGraph_slow(set<int> *previous_assign) {
    vector<node>::iterator nit, toit, fromit;
    vector<edge>::iterator eit;
    vector<int>::iterator eid_it;
    int nid, toid, fromid;
    int maxdiffpsi;
    vector<int> fassigned_edges;
    int remaining_assign=0, total_assign=0;

    initLists();

    for (nit=nodes.begin(), nid=0; nit!=nodes.end(), nid<noA; ++nit, ++nid) {
        maxdiffpsi=0;
        nit->best_el.idx=-1;
        getNextList(nid, nit->best_el);
        while (nit->best_el.idx!=-1 && nit->best_el.value<=nit->mindist) {
            toid=nit->best_el.idx;

            if (previous_assign[nid].find(toid)==previous_assign[nid].end()) { // not assigned
                toit=nodes.begin()+toid;
                insertEdge(nid, toid, nit->best_el.value, nit, toit);

                eit=edges.end()-1;
                if (calCost(eit, nit, toit)<0) {
                    if (maxdiffpsi<nit->psi-eit->weight)
                        maxdiffpsi=nit->psi-eit->weight;
                    //toit->psi=nit->psi-eit->weight;
                }
                //if (toit->isUpdatedfromPrevious) {
                //    toit->psi=nit->psi-eit->weight;
                //} else {
                //    nit->psi=toit->psi+eit->weight;
                //}
            }

            nit->best_el.idx=-1;
            getNextList(nid, nit->best_el);
        }

        nit->psi=nit->psi-maxdiffpsi;
        for (eid_it=nit->rE.begin(); eid_it!=nit->rE.end(); ++eid_it) {
            eit=edges.begin()+*eid_it;
            fromid=get_toid(nid, eit);
            fromit=nodes.begin()+fromid;
            if (calCost(eit, fromit, nit)<0)
                fassigned_edges.push_back(*eid_it);
        }
        total_assign+=nit->rE.size();

        // reverse edge
        for (eid_it=fassigned_edges.begin(); eid_it!=fassigned_edges.end(); ++eid_it)
            removeAssign(*eid_it, nid);

        fassigned_edges.clear();

        remaining_assign+=nit->rE.size();

        previous_assign[nid].clear();
    }

    printf("%d/%d\n", remaining_assign, total_assign);

#ifdef DEBUG
    for (nit=nodes.begin(), nid=0; nit!=nodes.end(); ++nit, ++nid) {
        for (eid_it=nit->E.begin(); eid_it!=nit->E.end(); ++eid_it) {
            eit=edges.begin()+*eid_it;
            toid=get_toid(nid, eit);
            toit=nodes.begin()+toid;

            if (calCost(eit, nit, toit)<0) {
                printf("Negative edges found (fromid=%d toid=%d eid=%d cost=%d)!\n", nid, toid, *eid_it, calCost(eit, nit, toit));
            }
        }
    }
#endif
}

void COptUpdate::updateSIA_slow() {
    // Stage 2
    runSIA();
}

void COptUpdate::reformSubGraph(set<int> *previous_assign) {
    vector<node>::iterator nit, toit, fromit;
    vector<edge>::iterator eit;
    vector<int>::iterator eid_it;
    int nid, toid, fromid;
    int maxdiffpsi;
    vector<int> fassigned_edges;

    initLists();

    for (nit=nodes.begin(), nid=0; nid<noA; ++nit, ++nid) {
        maxdiffpsi=0;
        nit->best_el.idx=-1;
        getNextList(nid, nit->best_el);
        while (nit->best_el.idx!=-1 && nit->best_el.value<=nit->mindist) {
            toid=nit->best_el.idx;

            if (previous_assign[nid].find(toid)==previous_assign[nid].end()) { // not assigned
                toit=nodes.begin()+toid;
                insertEdge(nid, toid, nit->best_el.value, nit, toit);

                eit=edges.end()-1;
                if (calCost(eit, nit, toit)<0) {
                    toit->psi=nit->psi-eit->weight;
                }
                //if (toit->isUpdatedfromPrevious) {
                //    toit->psi=nit->psi-eit->weight;
                //} else {
                //    nit->psi=toit->psi+eit->weight;
                //}
            }

            nit->best_el.idx=-1;
            getNextList(nid, nit->best_el);
        }

        previous_assign[nid].clear();
    }

#ifdef DEBUG
    for (nit=nodes.begin(), nid=0; nit!=nodes.end(); ++nit, ++nid) {
        for (eid_it=nit->E.begin(); eid_it!=nit->E.end(); ++eid_it) {
            eit=edges.begin()+*eid_it;
            toid=get_toid(nid, eit);
            toit=nodes.begin()+toid;

            if (calCost(eit, nit, toit)<0) {
                printf("Negative edges found (fromid=%d toid=%d eid=%d cost=%d)!\n", nid, toid, *eid_it, calCost(eit, nit, toit));
            }
        }
    }
#endif
}

/*
void COptUpdate::updateSIA() {
    // create artificial node
    node n;
    n.isUpdated=false;
    n.isUpdatedfromPrevious=false;
    n.psi=n.mindist=0;
    n.mineid=n.loop=-1;
    n.best_el.idx=-1;
    n.maxflow=Bcap*noB;

    nodes.push_back(n);
    vector<node>::iterator nit, toit;
    nit=nodes.begin()+(nodes.size()-1);

    for (int i=0; i<noB; ++i) {
        toit=nodes.begin()+i+noA;
        insertEdge(nodes.size()-1, i+noA, 0, nit, toit);
    }

    int count=0;
    for (int i=0; i<noA; ++i)
        if (nodes[i].flow>0)
            count++;

    int curdijkloop=0;
    // Stage 2
    while (!updates.empty()) {
        // maximize the capacity of the artificial node
        nodes.rbegin()->flow=Bcap*noB-Acap*noA;

        handleaffectedB(curdijkloop);

        updates.clear();

        fixNewPsi(curdijkloop);
    }

    count=0;
    for (int i=0; i<noA; ++i)
        if (nodes[i].flow>0)
            count++;

    // Stage 3
    runSIA(curdijkloop);
}
*/

void COptUpdate::handleaffectedB(int &curdijkloop) {
    vector<UpdateInfo>::iterator uit;
    vector<int> visited;
    int ret_aid;
    int sentflows=0;

    int i;
    int nid, fromid, tmpid, cost;
    vector<int>::iterator reid_it;
    vector<node>::iterator nit, fromit;
    vector<edge>::iterator eit;

    for (uit=updates.begin(); uit!=updates.end(); ++uit) {
        nit=nodes.begin()+uit->id;

        if (nodes[uit->id].flow==0 || nodes[uit->id].psi==0)
            //if (nodes[uit->id].flow==0)
            continue;

        initDijkstra(uit->id, curdijkloop, visited);
        ret_aid=compDijkstra_rE(curdijkloop, visited);

        if (ret_aid>-1) {
            // augment the flow
            int nopaths=augmentFlow_rE(ret_aid);

            // update potential
            updatePotential_rE(nodes[ret_aid].mindist, visited);

            // count the number of flows sent
            sentflows+=nodes[ret_aid].augflow;
        }
    }

    printf("rebuild psi\n");
    rebuildPsi(curdijkloop);
}

void COptUpdate::rebuildPsi(int &curdijkloop) {
    vector<node>::iterator nit, toit;
    vector<edge>::iterator eit;
    vector<int>::iterator eid_it;
    int nid, toid;
    queue<int> check;

    ++curdijkloop;

    // init A
    for (nit=nodes.begin(); nit!=nodes.begin()+noA; ++nit)
        nit->psi=0;

    // init B
    for (nit=nodes.begin()+noA, nid=noA; nit!=nodes.begin()+(noA+noB); ++nit, ++nid) {
        nit->psi=0;

        for (eid_it=nit->E.begin(); eid_it!=nit->E.end(); ++eid_it) {
            eit=edges.begin()+*eid_it;
            toid=get_toid(nid, eit);
            toit=nodes.begin()+toid;

            if (toid==noA+noB)
                continue;

            if (calCost(eit, nit, toit)<0) {
                toit->psi=eit->weight+nit->psi;
                if (toit->loop<curdijkloop || !toit->isCheck) {
                    toit->isCheck=true;
                    toit->loop=curdijkloop;
                    check.push(toid);
                }
            }
        }
    }

    while (!check.empty()) {
        nid=check.front();
        nit=nodes.begin()+nid;
        nit->isCheck=false;
        check.pop();

        for (eid_it=nit->E.begin(); eid_it!=nit->E.end(); ++eid_it) {
            eit=edges.begin()+*eid_it;
            toid=get_toid(nid, eit);
            toit=nodes.begin()+toid;

            if (toid==noA+noB)
                continue;

            if (calCost(eit, nit, toit)<0) {
                toit->psi=(toit->maxflow>0)? nit->psi+eit->weight : nit->psi-eit->weight;
                if (toit->loop<curdijkloop || !toit->isCheck) {
                    toit->isCheck=true;
                    toit->loop=curdijkloop;
                    check.push(toid);
                }
            }
        }
    }
}

void COptUpdate::fixNewPsi(int &curdijkloop) {
    vector<node>::iterator nit, toit, nit2;
    vector<edge>::iterator eit;
    vector<int>::iterator eid_it;
    int nid, toid;
    UpdateInfo update;
    queue<int> remove_eid;

    ++curdijkloop;

    // init A
    for (nit=nodes.begin(), nid=0; nit!=nodes.begin()+noA; ++nit, ++nid) {
        while (nit->best_el.idx!=-1 && nit->best_el.value<nit->psi) {
            toid=nit->best_el.idx;
            toit=nodes.begin()+toid;
            insertEdge(nid, toid, nit->best_el.value, nit, toit);
            eit=edges.begin()+(edges.size()-1);

            if (calCost(eit, nit, toit)<0) {
                toit->psi=nit->psi-eit->weight;
                if (toit->loop<curdijkloop) {
                    update.id=toid;
                    updates.push_back(update);
                    toit->loop=curdijkloop;

                    // remove assign
                    for (eid_it=toit->E.begin(); eid_it!=toit->E.end(); ++eid_it)
                        remove_eid.push(*eid_it);

                    while (!remove_eid.empty()) {
                        eit=edges.begin()+remove_eid.front();
                        nit2=nodes.begin()+get_toid(toid, eit);
                        if (toit->psi<nit2->psi-eit->weight)
                            toit->psi=nit2->psi-eit->weight;

                        removeAssign(remove_eid.front(), get_toid(toid, eit));
                        remove_eid.pop();
                    }
                }
            }

            nit->best_el.idx=-1;
            getNextList(nid, nit->best_el);
        }
    }
}

int COptUpdate::compDijkstra_rE(int &curloop, vector<int> &visited) {
    vector<int>::iterator eid_it;
    vector<node>::iterator fromit;
    vector<node>::iterator toit;
    vector<edge>::iterator eit;
    int curid, toid;
    int mindist, cost, ret, augflow, tmp;

    while (dijkH.size()>0) {
        dijkH.dequeue(curid, mindist);
        visited.push_back(curid);

        fromit=nodes.begin()+curid;

        // reach destination
        if (fromit->flow>0) return curid;

        for (eid_it=fromit->rE.begin(); eid_it!=fromit->rE.end(); ++eid_it) {
            eit=edges.begin()+*eid_it;
            toid=get_toid(curid, eit);
            toit=nodes.begin()+toid;

            if ((augflow=getNoAugment(eit, toit, fromit))>-1) {
#ifdef DEBUG
                if (calCost(eit, toit, fromit)<0) {
                    printf("bug %d %d %d\n", curid, toid, *eid_it);
                }
#endif

                if ((ret=updateMinDist(curloop, *eid_it, eit, fromit, toit, augflow, calCost(eit, toit, fromit)))==NEWNODE) {
                    dijkH.enqueue(toid, toit->mindist);
                } else if (ret==UPDATE) {
                    dijkH.updatequeue(toid, toit->mindist);
                }
            }
        }
    }

    return -1;
}

int COptUpdate::augmentFlow_rE(int lastid) {
    vector<edge>::iterator eit;
    vector<node>::iterator nit=nodes.begin()+lastid;
    vector<node>::iterator toit;
    int augflow=nit->augflow;
    int nopaths=0;

    nit->flow=nit->flow-augflow;

    while (nit->mineid!=-1) {
        toit=nit;
        eit=edges.begin()+nit->mineid;
        eit->flow=((nit->maxflow>0)? eit->flow+augflow : eit->flow-augflow);

        lastid=get_toid(lastid, eit);
        nit=nodes.begin()+lastid;

        int from_eidx=eit->node_eidx;
        int to_reidx=eit->node_reidx;

        insert_eid(toit->rE, eit->node_reidx, toit->mineid);
        remove_eid(nit->rE, to_reidx, edges[*toit->rE.rbegin()].node_reidx);

        insert_eid(nit->E, eit->node_eidx, toit->mineid);
        remove_eid(toit->E, from_eidx, edges[*toit->E.rbegin()].node_eidx);

        ++nopaths;
    }

    nit=nodes.begin()+lastid;
    nit->flow=nit->flow+augflow;

    return nopaths;
}


void COptUpdate::updatePotential_rE(int mindist, vector<int> &visited) {
    vector<int>::iterator it;
    vector<node>::iterator nit;
    for (it=visited.begin(); it!=visited.end(); ++it) {
        nit=nodes.begin()+*it;
        if (nit->mindist<mindist)
            nit->psi=nit->psi+nit->mindist-mindist;
    }
}
