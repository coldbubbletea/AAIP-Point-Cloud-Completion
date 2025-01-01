#include "OptAssign.h"
#include "SpatialCCA.h"
#include "utilities.h"

#include "SpatialCCA.h"
#include "rtree/rtree.h"
#include "rtree/blk_file.h"
#include<string>
int capQ[10005];
int capP[10005];

char pointfile[255];


RTree* buildRtree()
{

	ifstream real_point("../data/gt 6.txt");
	string pointfile_s = "../data/pointfile2";
	strcpy(pointfile, pointfile_s.c_str());
	Cache* c = new Cache(20000, 1024);
	RTree* a = new RTree(pointfile, 1024, c, 3);
	for (auto i = 1; i <= 8000; ++i)
	{
		Entry* b = new Entry(3, nullptr);
		float c1, c2, c3;
		real_point >> c1;
		//cout << "c1: "<<c1 << endl;
		real_point >> c2;
		real_point >> c3;
		//cout << "c1: " << c2 << endl;
		//cout << "c1: " << c3 << endl;
		b->num_data = 1;
		b->son = i;
		b->bounces[0] = b->bounces[1] = c1;
		b->bounces[2] = b->bounces[3] = c2;
		b->bounces[4] = b->bounces[5] = c3;
		
		a->insert(b);
	}
	real_point.close();
	return a;

}

void main(int argc, char* argv[]) {
    
    CSpatialCCA test;

    const char *qrytype, *objtype, *processtype;
    int qrysize, objsize, qryk, objk;
    int output;
    int no_qryupdates, no_objupdates;
    float randomized_ratio, update_ratio, insert_ratio;

    ConfigType cr;
    AddConfigFromFile(cr,"config.ini");
    AddConfigFromCmdLine(cr,argc,argv);

    output=getConfigInt("output", cr);
    qrytype=getConfigStr("qrytype", cr);
    objtype=getConfigStr("objtype", cr);
    qrysize=getConfigInt("qrysize", cr);
    objsize=getConfigInt("objsize", cr);
    qryk=getConfigInt("qryk", cr);
    objk=getConfigInt("objk", cr);
    processtype=getConfigStr("process", cr);

	int dom = getConfigInt("dom", cr);
    
    no_qryupdates=getConfigInt("no_qryupdates", cr);
    no_objupdates=getConfigInt("no_objupdates", cr);
    randomized_ratio=getConfigFloat("randomized_ratio", cr);
    update_ratio=getConfigFloat("update_ratio", cr);
    insert_ratio=getConfigFloat("insert_ratio", cr);

    char qryfile[255], objfile[255], assignfile[255];
    double start_secs;
	RTree* rt= buildRtree();
    sprintf(qryfile, "../data/q%s_d2_q%d_dom%d", qrytype, qrysize, dom);
    sprintf(objfile, "../data/r%c2dc%c_d2_o%d_dom%d", objtype[0], objtype[1], objsize, dom);

    if (strcmp(processtype, "o")==0) {
        printf("%s\n%s\n", qryfile, objfile);

        test.readQry(qryfile, qrysize);
		cout << "happy" << endl;
		test.readObj(rt, objsize);
		
		
		cout << "pass pass" << endl;
        //test.readObj(objfile, objsize);
		FILE* capFile = fopen("../data/cap_temp.txt", "r");
		int ca = 0, cb = 0;
		for (int i = 0; i < 8000; i++)
		{
			//fscanf(capFile, "%d", &capQ[i]);
			ca += 1;
		}
		for (int i = 1; i <= 1; i++)
		{
			for (int j = 0; j < 8000; j++)
			{
				//fscanf(capFile, "%d", &capP[j]);
				//fscanf(capFile, "%d", &capP[j]);
				cb += 1;
			}
			cout << ca << "###" << cb << endl;
			start_secs=gettime();
			cout << "~" << endl;
			double t1 = clock();
			test.initSIA(SIA, qrysize, objsize, capQ, capP);
			cout << "pass pass 2" << endl;
			rt->load_root();
			test.runSIA();
			
			cout << clock() - t1;
			double res = IntToFloat(test.getAssignCost())/2000;
			printf("total cost: %lf\n", res);
			
			//test.printstat(processtype, objfile, output);
			system("pause");
			cout << "@" << endl;
		}
		fclose(capFile);
    }

	printf("total execution time: %lf\n", gettime()-start_secs);
}