#include "SpatialCCA.h"
#include "rtree/rtree.h"
#include<iostream>
#include<utility>
#include<vector>
#include<math.h>
#include<string>
#include<sstream>
#include <iostream>
#include <string>
#include <cstring>
/**
* Test the correctness of RTree 
* Verify the NNs result correctness by comparing 
* We test the top-10 NNs and if all pass then this query will labeled pass
* then we count the pass rate = pass_query / num_of_query
* use brute force as the base ans
* @ hermes.cheung 2022-12-12
* @ nid,nit is the index and iterator of query 
* @ nq is num of test query 
**/
void CSpatialCCA::testRTree()
{
    auto nit = nodes.begin();
    auto i = 0;
    queue<int> free;
    auto nq = 1000;
    auto pass_ans = 0;
    
    vector<int> bf;
    for (auto i = 0; i < noA; ++i, ++nit)
        if (nit->flow > 0)
            free.push(i);
    initLists();
    auto tag = 1;
    auto tt = 1;
    while (nq--)
    {
        
        cout << "nq: " << nq << endl;
        
        
        auto nid = free.front();
        cout << "nid" << nid << endl;
        free.pop();
        nit = nodes.begin() + nid;
        int cnt = 10;
        int top = 0;
        while (--cnt)
        {
            tt++;
            getNextList(nid, nit->best_el);
            // brute force vector
            if (nit->best_el.idx == bf_v[nid][top].first)
            {
                tag++;
               
                top++;
            }
            else
            {
                
                cout << "RT idx: " << nit->best_el.idx ; cout << "           brute force idx: " << bf_v[nid][top].first << endl;
                cout << "RT value: " << nit->best_el.value; cout << "        brute force value: " << bf_v[nid][top].second << endl;
                
                top++;
            }

                
        }
        if (tag) pass_ans++;
        cout << "tag: " << tag << endl;
    }
   
    cout << "pass_rate: " << double(tag * 1.0) / double(tt * 1.0)*100 <<"%" << endl;
   
}
CSpatialCCA::CSpatialCCA(void)
{
}

CSpatialCCA::~CSpatialCCA(void)
{
}

void CSpatialCCA::readQry( int qrysize,int  objsize,float* v) {
	int iDim;
	float tmp;
	Qry.resize(qrysize);
	QryH.resize(qrysize);
    iDim = 3;
	for (int i=0; i< qrysize; ++i) {
		for (int j=0; j<iDim; ++j) {
            //    cout<<v[i*3+0]<<endl;
                Qry[i].v[j] = v[i*3+j];

		}
	}
    


}





void CSpatialCCA::readObj(RTree* RT,char* objfile, int objsize) {
    string pre_file_title="RTreeData/RTreeData1024/";
    string file_str(objfile);
    pre_file_title+=file_str;
   
   // cout<<pre_file_title<<endl;
    char* c_str = new char[pre_file_title.length() + 1]; // Allocate memory for the character array
    std::strcpy(c_str, pre_file_title.c_str()); 
    
    rt = new RTree(c_str, 0);
	// Load rtree
	rt->load_root();
    //cout<<"yes"<<endl;
   /// cout<<endl<<"goto"<<endl;
	// initial vObjects
	//objsize=rt->num_of_data;
   // cout << "objsize: " << objsize << endl;
    delete[] c_str;
    RT=rt;
    //this->rt = rt;
}


void CSpatialCCA::readObjbounces(char* objfile)
{
    FILE *fd=fopen(objfile,"rb");

	Obj.resize(noB);

    for (int i=0; i<noB; ++i)
    {
        float fTemp;
        for (int j=0; j<3; ++j)
        {
            fTemp=0.0f;
            fread(&fTemp, sizeof(float), 1, fd);


            Obj[i].v[2*j]=Obj[i].v[2*j+1]=(float)fTemp;
        }

        fread(&fTemp, sizeof(float), 1, fd);
    }

    fclose(fd);
}

void CSpatialCCA::initLists() {
    //cout << "enter" << endl;
	vector<Point2D>::iterator qit;
    
	vector<Heap>::iterator qitH;
	int qid;
    int cnt = 0;
	for (qit=Qry.begin(), qitH=QryH.begin(), qid=0; qit!=Qry.end(), qitH!=QryH.end(); ++qit, ++qitH, ++qid) {
        
       // cout << "cnt: " << ++cnt << endl;
       // cout << "start of a query" << endl;
		while (qitH->size()>0) qitH->pop();	// clear the heap first
       // cout << "root: " << rt->root << endl;
		RTNode* cur_node=rt->root_ptr;
		for (int i=0;i<cur_node->num_entries;i++) {
            
            //cout<<"i"<<endl;
			HeapEntry he;
			he.level=cur_node->level;
			he.son=cur_node->entries[i].son;
			he.key=MINDIST_nosqrt(*qit,cur_node->entries[i].bounces);
			//memcpy(he.bounces, cur_node->entries[i].bounces, sizeof(float)*2*rt->dimension);
			qitH->push(he);
		}
	}
}

//bool CSpatialCCA::getNextList(int id, fHeap<int>::elem &el) { 
//	if (el.idx!=-1)
//		return true;
//
//	vector<Point2D>::iterator qit = Qry.begin()+id;
//	vector<Heap>::iterator qitH = QryH.begin()+id;
//
//	while (qitH->size()>0) {
//		HeapEntry he=qitH->top();	// dequeue next entry
//		qitH->pop();
//		if (he.level!= 0) 		{			
//			RTNode *rtn=new RTNode(rt,he.son);
//			for (int i=0;i<rtn->num_entries;i++) {
//				HeapEntry new_he;
//				new_he.key=MINDIST_nosqrt(*qit,rtn->entries[i].bounces);
//				new_he.level=rtn->level;
//				new_he.son=rtn->entries[i].son;
//				qitH->push(new_he);
//			}
//			delete rtn;
//			rtn=NULL;
//		} else {
//			el.idx=he.son+noA;
//			el.value=FloatToInt(he.key);
//			return true;
//		}
//	}
//
//	return false;
//}

bool CSpatialCCA::getNextList(int id, fHeap<int>::elem &el) { 
	/*if (el.idx!=-1)
		return true;*/

	vector<Point2D>::iterator qit = Qry.begin()+id;
	vector<Heap>::iterator qitH = QryH.begin()+id;

	while (qitH->size()>0) {
		HeapEntry he=qitH->top();	// dequeue next entry
		qitH->pop();
		if (he.level!= 0) 		{			
			RTNode *rtn=new RTNode(rt,he.son);
			for (int i=0;i<rtn->num_entries;i++) {
				HeapEntry new_he;
				new_he.key=MINDIST_nosqrt(*qit,rtn->entries[i].bounces);
				new_he.level=rtn->level;
				new_he.son=rtn->entries[i].son;
				qitH->push(new_he);
			}
			delete rtn;
			rtn=NULL;
		} else {
			map<int,int>::iterator mit;
			if ((mit=nodes[id].preAssign.find(he.son+noA))==nodes[id].preAssign.end()) {
				el.idx=he.son+noA;
				el.value=FloatToInt(he.key);
				return true;
			} else {
				edges[mit->second].weight=FloatToInt(he.key);
			}
		}
	}

	return false;
}

float CSpatialCCA::MINDIST_nosqrt(Point2D &point, float *bounces, int dim) {
    
	float r,sum=0.0;
	for (int i = 0; i < dim; i++) {
       // std::cout << "query: " << point.v[i] << endl;
        //cout << "tree1: " << bounces[2 * i + 1] << endl;
       // cout << "tree1: " << bounces[2 * i] << endl;
		r=max(bounces[2*i]-point.v[i], point.v[i]-bounces[2*i+1]);
		if (r>0) sum += r*r;
	}
   // cout << "sum: " << sum << endl;
   // cout << "end compare" << endl;

	return sum;
}

void CSpatialCCA::printstat(const char *processtype, char *objfile, int output) {
	vector<edge>::iterator eit;

	if (strcmp(processtype, "c")==0) {
		sprintf(objfile, "%s.txt", objfile);
		readObjbounces(objfile);
	}

    if (output>0) {
        vector<node>::iterator nit, toit;
        vector<int>::iterator eidit;
		int i, toid;

		if (strcmp(processtype, "i")!=0) {
			for (nit=nodes.begin()+noA, i=noA; i<noA+noB; ++nit, ++i) {
				for (eidit=nit->E.begin(); eidit!=nit->E.end(); ++eidit) {
					eit=edges.begin()+*eidit;
					toid=get_toid(i, eit);
					toit=nodes.begin()+toid;
					toit->rE.push_back(*eidit);

					if (strcmp(processtype, "c")==0) {
						if (eit->weight==0)
							eit->weight=FloatToInt(MINDIST_nosqrt(Qry[toid],Obj[i].v));
					}
				}
			}
		}
        for (nit=nodes.begin(), i=0; i<Qry.size(); ++nit, ++i) {
            printf("%d: ", i);
            for (eidit=nit->rE.begin(); eidit!=nit->rE.end(); ++eidit) {
                eit=edges.begin()+*eidit;
                if (eit->flow>0)
                    printf("%d flow to %d (w: %d) ", eit->flow, eit->toid-noA, eit->weight);
            }
            printf("\n");
        }
    }

	printf("\nno of edges: %d\n", edges.size());
	//printf("total assignment cost: %lf\n", IntToFloat(getAssignCost()));
}


void CSpatialCCA::readObjUpdates(char* objupdatefile) {
    sprintf(objupdatefile, "%s.txt", objupdatefile);
    FILE *fd=fopen(objupdatefile,"r");
   
    if (fd==NULL) { printf("No updates file exists!\n"); exit(-1); }

    int i;
    UpdateInfo uinfo;
    char update, type;
    int transform;

    fscanf(fd, "%d\n", &no_objupdates);
    for (i=0; i<no_objupdates; ++i) {
        fscanf(fd, "%c %c %d", &update, &type, &uinfo.id);
        
        // convert object id to node id
        uinfo.id=uinfo.id+Qry.size();
        
        switch (type) {
            case 'q': uinfo.type=A; break;
            case 'o': uinfo.type=B; break;
        }

        switch (update) {
            case 'u': uinfo.update=UP; break;
            case 'i': uinfo.update=INS; break;
            case 'd': uinfo.update=DEL; break;
        }

        if (update=='d') {
            fscanf(fd, "%d\n", &transform);
            delid[transform+Qry.size()]=uinfo.id;
        } else {
            fscanf(fd, "\n");
        }

        nodes[uinfo.id].isUpdatedfromPrevious=true;

        if (uinfo.type==B && uinfo.update!=DEL)
            updates.push_back(uinfo);
    }
    fclose(fd);
}

void CSpatialCCA::readAssignment(char* assignfile) {
    FILE *fd=fopen(assignfile,"rb");

    if (fd==NULL) { printf("No assignment file exists!\n"); exit(-1); }

    vector<node>::iterator nit, nit2, toit;
    map<int,int>::iterator mit;
    int nid, k, dist, toid;
    int maxpsi=0, psi;
    int qryid, objid;
    int num_assign;

    for (nit=nodes.begin(), nid=0; nit!=nodes.end(); ++nit, ++nid) {
        fread(&psi, sizeof(int), 1, fd);

        if ((mit=delid.find(nid))!=delid.end())
            nit2=nodes.begin()+mit->second;
        else
            nit2=nit;

        nit2->psi=psi;

        if (nid<Qry.size() && nit->psi>maxpsi)
            maxpsi=psi;
    }

    // [IMPORTANT] init previous_assign
    set<int>* previous_assign=new set<int>[Qry.size()];

    while (!feof(fd)) {
        fread(&qryid, sizeof(int), 1, fd);
        if (qryid==-1)
            break;

        nid=qryid;

        nit=nodes.begin()+qryid;
        fread(&nit->mindist, sizeof(int), 1, fd);
        fread(&num_assign, sizeof(int), 1, fd);

        // [IMPORTANT] set the largest psi to mindist
        nit->mindist=nit->psi;

        for (int i=0; i<num_assign; ++i) {
            fread(&objid, sizeof(int), 1, fd);
            fread(&k, sizeof(int), 1, fd);
            fread(&dist, sizeof(int), 1, fd);
            
            if ((mit=delid.find(objid))!=delid.end())
                toid=mit->second+Qry.size();
            else
                toid=(objid+Qry.size());

            toit=nodes.begin()+toid;

            if (!nit->isUpdatedfromPrevious && !toit->isUpdatedfromPrevious) {
				insertAssign(nid, toid, dist, nit, toit);
                //insertEdge(toid, nid, dist, toit, nit);
                //--nit->flow;
                //++toit->flow;

                // [IMPORTANT] insert toid into previous assign
                previous_assign[nid].insert(toid);
			}
        }
    }

    fclose(fd);

    reformSubGraph(previous_assign);

    delete[] previous_assign;
}
