#include "OptAssign.h"
#include <map>
#include <chrono>
#include <set>
#include <utility>
#include <limits.h>
#include<chrono>
#include<iostream>
#include<fstream>
COptAssign::~COptAssign(void)
{
	vector<node>().swap(nodes);
	vector<edge>().swap(edges);
}

void COptAssign::initSIA(int runtype, int noA, int noB, int Acap[], int Bcap[], int edgecap)
{
	// nodes initialization
	this->noA = noA;
	this->noB = noB;
	// this->Acap=Acap;
	// this->Bcap=Bcap;
	this->edgecap = edgecap;
	this->runtype = runtype;

	int sumcapA = 0, sumcapB = 0;
	for (int i = 0; i < noA; i++)
	{
		sumcapA += Acap[i];
	}
	for (int i = 0; i < noB; i++)
	{
		sumcapB += Bcap[i];
	}

	totalflows = min(sumcapA, sumcapB);
	nodes.clear();
	edges.clear();
	updateH.clear();
	dijkH.clear();
	maxpsi = 0; // ?
	nodes.reserve(noA + noB);

	node n;
	n.isUpdated = false;
	n.isUpdatedfromPrevious = false;
	n.psi = n.mindist = 0;
	n.mineid = n.loop = -1;
	n.best_el.idx = -1;

	for (int i = 0; i < noA; ++i)
	{
		//cout<<"Acap[i]: "<<Acap[i]<<endl;
		n.flow = n.maxflow = Acap[i];
		nodes.push_back(n);
	}

	for (int i = 0; i < noB; ++i)
	{
		n.flow = n.maxflow = -Bcap[i];
		nodes.push_back(n);
	}
		// Lp-norm optimization heuristics
	// for (int i = 0; i < noA; i++)
	// {
	// 	int tflow = std::min(nodes[i].flow, -nodes[i+noA].flow);
	// 	nodes[i].flow -= tflow;
	// 	nodes[i+noA].flow += tflow;
	// 	totalflows -= tflow;
	// }

	assignEdges.clear(); // for upper bound maintaineance- Yu
	asgn.resize(noA + noB);
	for (int i = 0; i < noA + noB; i++)
	{
		asgn[i].clear();
	}
}
	

int COptAssign::removeFeasibleAssign(int a, int b, int delta)
{
	// a - augment path start, b - augment path end
	// free delta capacity from both a and b
	// return reduced cost from the feasible solution that we are maintaining

	int removedCost = 0;
	for ( ; asgnHead[b] < asgn[b].size(); asgnHead[b]++)
	{
		if (delta == 0) break;
		int pid = *(asgn[b].begin() + asgnHead[b]);
		int ex = std::min(delta, assignEdges[pid].cap);
		if (ex == 0)
			continue;
		int tmpex = ex;
		int y = assignEdges[pid].from;
		for ( ; asgnHead[a] < asgn[a].size(); asgnHead[a]++)
		{
			if (tmpex == 0) break;
			int qid = *(asgn[a].begin() + asgnHead[a]);
			if (pid == qid)
			{
				continue;
			}
			int tmp = std::min(tmpex, assignEdges[qid].cap);
			if (tmp == 0)
				continue;
			tmpex -= tmp;
			assignEdges[qid].cap -= tmp;
			removedCost += assignEdges[qid].cost * tmp;
			
			int x = assignEdges[qid].to;
			//int tcost = TIMESFACTOR * costs[y][x - noA];
			int tcost = costs[y][x - noA];
			assignEdges.push_back(AssignE(y, x, tcost, tmp));
			int eid = assignEdges.size() - 1;
			asgn[y].push_back(eid);
			asgn[x].push_back(eid);
			removedCost -= tcost * tmp;
		}
		delta -= ex;
		assignEdges[pid].cap -= ex;
		removedCost += assignEdges[pid].cost * ex;
	}
	return removedCost;
}
int COptAssign::runSIA(int * is_conflict,int * curloopsent , int curdijkloop ) {
	
	auto test_s = std::chrono::system_clock::now();
	int i, j;
	vector<node>::iterator nit;
	queue<int> free;
	int retid;
	int tmp;
	int path_cost;
	fHeap<int> globalH;
	vector<int> visited;
	vector<int>::iterator it;
	int availablePath;
	int sentflows=0;
	int curloop=0;
	vector<int> updated_nodes;
	long long current_running_time_milliseconds=0;

	// free Nodes
	// for (i=0, nit=nodes.begin(); i<noA; ++i, ++nit)
	// 	if (nit->flow>0)
	// 		free.push(i);

	

	
	/******************** For naive SSPA ********************/
	//for (int nid=0; nid<noA; ++nid) {
	//	vector<node>::iterator nit=nodes.begin()+nid;
	//	while (getNextList(nid, nit->best_el)) {
	//		vector<node>::iterator toit=nodes.begin()+nit->best_el.idx;
	//		insertEdge(nid, nit->best_el.idx, nit->best_el.value, nit, toit);
	//		nit->best_el.idx=-1;
	//	}
	//}
	/******************** For naive SSPA ********************/


	//feasiblecost = getFeasibleSolution(); // for upper bound
	//printf("#%f\n", IntToFloat(feasiblecost));
   
	currentcost = 0;
	remaincost = 0;
	currentflow = 0;
	int upperEst;
	std::map<int,std::set<int>> cflt_mp;
	int is_InSert[item_size];
	for(auto i=0;i<noA;++i) is_InSert[i]=0;
    //for(auto i=0;i<noA;++i) is_InSert[i]=0;
	// for(auto i=0;i<noA;++i)
	// {
	// 	//cout<<"conflict :"<<is_conflict[i]<<"\n"<<endl;
	// 	cflt_mp[is_conflict[i]].insert(i);
	// }
	// for(auto i=0;i<noA;++i)
	// {
	// 	if(cflt_mp[i].size()>1)
	// 	{
	// 		for(auto j=cflt_mp[i].begin();j!=cflt_mp[i].end();++j)
	// 		{
				
	// 			free.push(*j);
	// 			is_InSert[*j]=1;
	// 		}
			
	// 	}
	// }
	// int conflict_node_size=0;
    // conflict_node_size=free.size();
	// cout<<"conflict_node_size: "<<conflict_node_size<<endl;
	for(auto i=0;i<noA;++i)
	{
		if(!is_InSert[i]) free.push(i); 
	}
	//cout<<"free size"<<free.size()<<endl;
	if (runtype==SIA) initLists();
	vector<int> mindist(noA, 0);
	while (!free.empty() && sentflows<totalflows) 
	{
		int nid=free.front();
		auto start = std::chrono::system_clock::now();
		if (current_running_time_milliseconds >= shutdown_time_milliseconds)
		{
			// printf("\nmatching rate>70%\n");

			cout << "curloop"<< " " << current_running_time_milliseconds << endl;
			*curloopsent=curloop;
			break;
		}	
			free.pop();

		nit=nodes.begin()+nid;

		globalH.reset();
		if (getNextList(nid, nit->best_el))
			globalH.enqueue(nid, nit->best_el.value-nit->psi);
		maxpsi=nit->psi;

		//SIA_opt1(nid, nit);
		initDijkstra(nid, curdijkloop, visited);
		retid=compDijkstra_INC(curdijkloop, globalH, visited);
		while (retid==-1) {// || nodes[retid].mindist>globalH.getTopValue()-maxpsi) {
			availablePath=0;
			//while (availablePath<1)
			while (dijkH.size()==0 || dijkH.getTopValue()>globalH.getTopValue()) // Kamiru
				availablePath+=insertEdgefromHeap(curdijkloop, globalH);

			retid=updateDijkstra(curdijkloop, globalH, visited);
		}

		// augment the flow


		//printf("[%d %d] ", nid, retid);
		// update potential

		//est += (nit->maxflow - (nit->flow - augflows))  * (cur_path_cost - cur_dist[nid]);

		// --- augment a path from nid to retid

		int augflows=augmentFlow(retid, path_cost);

		sentflows+=augflows;
		currentflow = sentflows;
		
		vector<node>::iterator tit = nodes.begin() + retid;

		int cur_path_cost = path_cost + nit->psi - tit->psi; // from nit to tit

		int add_cost = augflows * cur_path_cost;
		
		//cout<<augflows<<" "<<cur_path_cost<<endl;
		
		currentcost += add_cost;
		
		remaincost += (cur_path_cost - mindist[nid]) * (nit->flow + augflows) - add_cost;
		mindist[nid] = cur_path_cost;

		int est = currentcost + remaincost; // quick bound

		upperEst = currentcost  + feasiblecost;

		
		if (IntToFloat(est) > 1e20) //Modified by Edison 
		{
			return est;
			//return -1; // prune
		}

		updatePotential(nodes[retid].mindist, visited);

		if (nodes[nid].flow>0) free.push(nid);

		++curloop;
		auto end = std::chrono::system_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		current_running_time_milliseconds += static_cast<long long int>(elapsed.count());
		
	}
	return curloop;
}

int COptAssign::insertEdgefromHeap(int curloop, fHeap<int> &globalH)
{
	int nid, toid, value, tmp;
	globalH.dequeue(nid, value);

	vector<node>::iterator nit = nodes.begin() + nid;
	vector<node>::iterator toit = nodes.begin() + nit->best_el.idx;
	vector<edge>::iterator eit;

	insertEdge(nid, nit->best_el.idx, nit->best_el.value, nit, toit);
	//   nit->E.push_back(edges.size());
	// toit->rE.push_back(edges.size());

	//   edge e;
	//   e.fromid=nid;
	//   e.toid=nit->best_el.idx;
	//   e.weight=nit->best_el.value;
	//   e.flow=0;
	// e.node_eidx=nit->E.size()-1;
	// e.node_reidx=toit->rE.size()-1;
	//   edges.push_back(e);

	// update dijkstra
	if (!toit->isUpdated)
	{
		// node[nid] has been popped already
		if (nit->loop == curloop && !dijkH.isExisted(nid, tmp))
		{
			eit = edges.end() - 1;
			int ret = updateMinDist(curloop, edges.size() - 1, eit, nit, toit, getNoAugment(edges.end() - 1, nit, toit), calCost(eit, nit, toit));
			updateHeaps(curloop, nit->best_el.idx, toit, ret);
		}
	}

	// clear the current best
	nit->best_el.idx = -1;
	if (getNextList(nid, nit->best_el))
	{
		globalH.enqueue(nid, (nit->mindist) + (nit->best_el.value) - (nit->psi));
		// cout << "nid: " << nid << endl;
		// cout << "nit->best_el.idx" << nit->best_el.idx<< endl;
	}

	return (toit->flow < 0) ? 1 : 0; //  if(flow<0) return 1 else return 0
}
void COptAssign::updateHeaps(int curloop, int toid, vector<node>::iterator toit, int ret)
{
	int tmp;
	if (ret == NEWNODE)
		// node[toid] has not been seen
		dijkH.enqueue(toid, toit->mindist);
	else if (ret == UPDATE)
	{
		// node[toid] has been seen
		if (dijkH.isExisted(toid, tmp))
			// node[toid] is in dijkH
			dijkH.updatequeue(toid, toit->mindist);
		else
		{
			// node[toid] has been popped
			if (!updateH.isExisted(toid, tmp))
				updateH.enqueue(toid, toit->mindist);
			else
				updateH.updatequeue(toid, toit->mindist);
			toit->isUpdated = true;
		}
	}
}

void COptAssign::initDijkstra(int nid, int &curloop, vector<int> &visited)
{
	++curloop;

	dijkH.reset2();
	dijkH.enqueue(nid, 0);
	visited.clear();

	vector<node>::iterator nit = nodes.begin() + nid;
	nit->mineid = -1;
	nit->mindist = 0;
	nit->loop = curloop;
	nit->augflow = abs(nit->flow);
}

int COptAssign::compDijkstra_INC(int &curloop, fHeap<int> &globalH, vector<int> &visited)
{
	vector<int>::iterator eid_it;
	vector<node>::iterator fromit;
	vector<node>::iterator toit;
	vector<edge>::iterator eit;
	int curid, toid;
	int mindist, cost, ret, augflow, tmp;

	while (dijkH.size()>0 && (globalH.size()==0 || dijkH.getTopValue()<=globalH.getTopValue()))  // Kamiru
	{ // Kamiru
		dijkH.dequeue(curid, mindist);
		visited.push_back(curid);

		fromit = nodes.begin() + curid;

		if (fromit->maxflow > 0 && !globalH.isExisted(curid, tmp))
		{
			if (getNextList(curid, fromit->best_el))
				globalH.enqueue(curid, fromit->mindist + fromit->best_el.value - fromit->psi);

			maxpsi = max(maxpsi, fromit->psi);
		}

		// reach destination
		if (fromit->flow < 0)
			return curid;

		for (eid_it = fromit->E.begin(); eid_it != fromit->E.end(); ++eid_it)
		{
			eit = edges.begin() + *eid_it;
			toid = get_toid(curid, eit);
			toit = nodes.begin() + toid;

			if (toid == noA + noB)
				continue;

			if ((augflow = getNoAugment(eit, fromit, toit)) > -1)
			{
#ifdef DEBUG
				if (calCost(eit, fromit, toit) < 0)
				{
					printf("bug %d %d %d\n", curid, toid, *eid_it);
				}
#endif

				if ((ret = updateMinDist(curloop, *eid_it, eit, fromit, toit, augflow, calCost(eit, fromit, toit))) == NEWNODE)
				{
					dijkH.enqueue(toid, toit->mindist);
				}
				else if (ret == UPDATE)
				{
					dijkH.updatequeue(toid, toit->mindist);
				}
			}
		}
	}

	return -1;
}

int COptAssign::updateDijkstra(int &curloop, fHeap<int> &globalH, vector<int> &visited)
{
	vector<int>::iterator eid_it;
	vector<edge>::iterator eit;
	vector<node>::iterator fromit, toit;
	int curid, toid, mindist, augflow, ret;

	while (updateH.size() > 0)
	{
		updateH.dequeue(curid, mindist);
		fromit = nodes.begin() + curid;
		fromit->isUpdated = false;

		for (eid_it = fromit->E.begin(); eid_it != fromit->E.end(); ++eid_it)
		{
			eit = edges.begin() + *eid_it;
			toid = get_toid(curid, eit);
			toit = nodes.begin() + toid;

			if ((augflow = getNoAugment(eit, fromit, toit)) > -1)
			{
				ret = updateMinDist(curloop, *eid_it, eit, fromit, toit, augflow, calCost(eit, fromit, toit));
				updateHeaps(curloop, toid, toit, ret);
			}
		}
	}

	return compDijkstra_INC(curloop, globalH, visited);
}

int COptAssign::updateMinDist(int curloop, int eid, vector<edge>::iterator eit, vector<node>::iterator fromit, vector<node>::iterator toit, int augflow, int cost)
{
	// int cost=calCost(eit, fromit, toit);
	if (toit->loop < curloop)
	{
		toit->mindist = fromit->mindist + cost;
		toit->mineid = eid;
		toit->augflow = min(fromit->augflow, augflow);
		toit->loop = curloop;
		return NEWNODE;
	}
	else if (toit->mindist > fromit->mindist + cost)
	{
		toit->mindist = fromit->mindist + cost;
		toit->mineid = eid;
		toit->augflow = min(fromit->augflow, augflow);
		return UPDATE;
	}
	else
		return NOUPDATE;
}

int COptAssign::augmentFlow(int lastid, int & cost) {
	vector<edge>::iterator eit;
	vector<node>::iterator nit;
	vector<node>::iterator toit;
	vector<int>::iterator it;
	int tempid;
	int max_augflow=INT_MAX;
	int nopaths=0;
	int edge_cap;

	cost = 0; // Yu

	// first phase: calculate the maximum flows being sent
	tempid=lastid;
	nit=nodes.begin()+tempid;
	if (max_augflow>abs(nit->flow))
		max_augflow=abs(nit->flow);
	while (nit->mineid!=-1) {
		toit=nit;
		eit=edges.begin()+nit->mineid;
		edge_cap=(nit->maxflow>0)? eit->flow : eit->maxflow-eit->flow;
		if (max_augflow>edge_cap)
			max_augflow=edge_cap;

		tempid=get_toid(tempid, eit);
		nit=nodes.begin()+tempid;
		cost += calCost(eit, nit, toit); // Yu
	}

	if (max_augflow>abs(nit->flow))
		max_augflow=abs(nit->flow);

	// second phase: update the flow graph
	tempid=lastid;
	nit=nodes.begin()+tempid;
	nit->flow=nit->flow+max_augflow;
	while (nit->mineid!=-1) {
		toit=nit;
		eit=edges.begin()+nit->mineid;
		eit->flow=(nit->maxflow>0)? eit->flow-max_augflow : eit->flow+max_augflow;

		// Yu
		

		tempid=get_toid(tempid, eit);
		nit=nodes.begin()+tempid;

		int from_eidx=eit->node_eidx;
		int to_reidx=eit->node_reidx;

		// insert an edge IFF it does not exist in E
		for (it=toit->E.begin(); it!=toit->E.end(); ++it)
			if (*it==toit->mineid)
				break;
		if (it==toit->E.end())
			insert_eid(toit->E, eit->node_eidx, toit->mineid);
		
		// remove an edge IFF there is no remaining capacity on the edge
		if ((nit->maxflow>0 && eit->flow==eit->maxflow) || // from A to B
			(nit->maxflow<0 && eit->flow==0))              // from B to A
			remove_eid(nit->E, from_eidx, edges[*nit->E.rbegin()].node_eidx);

		if (runtype==SIAUPDATE) {
			insert_eid(nit->rE, eit->node_reidx, toit->mineid);
			remove_eid(toit->rE, to_reidx, edges[*toit->rE.rbegin()].node_reidx);
		}

		++nopaths;
	}

	nit=nodes.begin()+tempid;
	nit->flow=nit->flow-max_augflow;

	return max_augflow;
}

void COptAssign::insertEdge(int fromid, int toid, int dist, vector<node>::iterator fromit, vector<node>::iterator toit)
{
	fromit->E.push_back(edges.size());
	if (runtype == SIAUPDATE)
		toit->rE.push_back(edges.size());

	edge e;
	if (fromid < toid)
	{
		e.fromid = fromid;
		e.toid = toid;
	}
	else
	{
		e.fromid = toid;
		e.toid = fromid;
	}
	e.weight = dist;
	e.flow = 0;
	e.node_eidx = fromit->E.size() - 1;
	e.node_reidx = toit->rE.size() - 1;
	edges.push_back(e);
}

void COptAssign::insert_eid(vector<int> &E, int &node_eidx, int eid)
{
	E.push_back(eid);
	node_eidx = E.size() - 1;
}

void COptAssign::remove_eid(vector<int> &E, int del_eidx, int &updated_eidx)
{
	if (E.size() - 1 > del_eidx)
	{
		E[del_eidx] = *E.rbegin();
		updated_eidx = del_eidx;
	}
	E.pop_back();
}

void COptAssign::updatePotential(int mindist, vector<int> &visited)
{
	vector<int>::iterator it;
	vector<node>::iterator nit;
	for (it = visited.begin(); it != visited.end(); ++it)
	{
		nit = nodes.begin() + *it;
		if (nit->mindist < mindist)
			nit->psi = nit->psi - nit->mindist + mindist;
	}
}

void COptAssign::SIA_opt1(int fromid, vector<node>::iterator fromit)
{
	vector<int>::iterator eid_it;
	vector<node>::iterator toit;
	vector<edge>::iterator eit;
	int aug;

	// getBestEdge(fromid, fromit, eit);

	for (eid_it = fromit->E.begin(); eid_it != fromit->E.end(); ++eid_it)
	{
		eit = edges.begin() + *eid_it;
		toit = nodes.begin() + get_toid(fromid, eit);

		// if ((aug=getNoAugment(fromid, eit, fromit, toit))>-1) {
		// }
		// if (toit->
	}
}

// inline int COptAssign::get_toid(int nid, vector<edge>::iterator eit) {
// 	return (eit->fromid==nid)? eit->toid : eit->fromid;
// };

// inline int COptAssign::getNoAugment(vector<edge>::iterator eit, vector<node>::iterator fromit, vector<node>::iterator toit) {
// 	// int flow=min(edgecap, min(fromit->maxflow,-toit->maxflow))-eit->flow; // THE bug!!!
// 	int flow=min(fromit->maxflow,-toit->maxflow)-eit->flow;
// 	if (fromit->maxflow>0 && flow>0)
// 		return flow;
// 	else if (fromit->maxflow<0 && eit->flow>0)
// 		return eit->flow;
// 	else return -1;
// };

// inline int COptAssign::calCost(vector<edge>::iterator eit, vector<node>::iterator fromit, vector<node>::iterator toit) {
// 	return (fromit->maxflow>0)? eit->weight-fromit->psi+toit->psi : -eit->weight-fromit->psi+toit->psi ;
// }

void COptAssign::insertAssign(int fromid, int toid, int dist, vector<node>::iterator fromit, vector<node>::iterator toit)
{
	insertEdge(toid, fromid, dist, toit, fromit);

	--fromit->flow;
	++toit->flow;

	++edges.rbegin()->flow;
}

void COptAssign::removeAssign(int eid, int toid)
{
	vector<edge>::iterator eit = edges.begin() + eid;
	vector<node>::iterator fromit = nodes.begin() + get_toid(toid, eit);
	vector<node>::iterator toit = nodes.begin() + toid;

	int from_eidx = eit->node_eidx;
	int to_reidx = eit->node_reidx;

	insert_eid(toit->E, eit->node_eidx, eid);
	remove_eid(fromit->E, from_eidx, edges[*fromit->E.rbegin()].node_eidx);

	insert_eid(fromit->rE, eit->node_reidx, eid);
	remove_eid(toit->rE, to_reidx, edges[*toit->rE.rbegin()].node_reidx);

	toit->flow = toit->flow + 1;
	fromit->flow = fromit->flow - 1;

	--eit->flow;
}

long long COptAssign::IPA(int *assignment, int *assignment_inv, float *price, int qrysize)
{
	int sum_non_zero_price=0;
	std::map<int,std::set<std::pair<int,float>>> mp;
	vector<edge>::iterator eit;
	float Assignment[item_size];
	#pragma omp parallel for
	for (auto i = 0; i < item_size; i++) Assignment[i] = 0;
		
	long long totalcost = 0;
    auto start = std::chrono::system_clock::now();
	for (eit = edges.begin(); eit != edges.end(); eit++)
	{
		// cout << "from: " << eit->fromid << " " << "to: " << eit->toid<<"weight: "<< eit->weight/5931.0f/100 <<"flow: "<<eit->flow << endl;

		if (eit->flow != 0)
			Assignment[eit->fromid]=eit->weight;  //record the d_kl
		else
			mp[eit->toid-item_size].insert(make_pair(eit->fromid,eit->weight)); //record s_k and dkj

	}
	for (auto i = 0; i < item_size; ++i)
	{
		
		
		for (auto j=mp[i].begin();j!=mp[i].end();j++)
		{
			
			int d_kj = j->second;
			int s_k = j->first;
			float d_kl = Assignment[s_k];
			float m=d_kl-d_kj;
			m=IntToFloat(m);
			price[i]=max(price[i],m);
			
		}
	}
	auto end = std::chrono::system_clock::now();
	
    #pragma omp parallel for
	for (auto i=0;i<item_size;++i) if(price[i]) sum_non_zero_price++;
	
	return 0;

	
}

void COptAssign::computeConstraint()
{
	vector<node>::iterator fromit, toit;
	int fromid, toid;

	toid = noA;
	toit = nodes.begin() + toid;

	/*
	for (int i=0; i<Bcap; ++i) {
		fromit=nodes.begin()+i;

		//insertAssign(i, toid, 0, fromit, toit);

		//fromit->preAssign.insert(pair<int,int>(toid, edges.size()-1));
	}
	*/

	runSIA();
}