#include <torch/extension.h>
#include <vector>
#include "OptAssign.h"
#include "SpatialCCA.h"
#include "utilities.h"
#include <chrono>
#include "SpatialCCA.h"
#include "rtree/rtree.h"
#include "rtree/blk_file.h"
#include <string>
#include <sstream>
#include <iostream>
#include <torch/extension.h>
#include <vector>
#include <iostream>
#include <iostream>
#include <fstream>
#include <cstdio>
using namespace chrono;
int capQ[100005];
int capP[100005];
int emd_cuda_forward(at::Tensor xyz1, at::Tensor xyz2, at::Tensor dist, at::Tensor assignment, at::Tensor price,
					 at::Tensor assignment_inv, at::Tensor bid, at::Tensor bid_increments, at::Tensor max_increments,
					 at::Tensor unass_idx, at::Tensor unass_cnt, at::Tensor unass_cnt_sum, at::Tensor cnt_tmp, at::Tensor max_idx, float eps, int iters);

int emd_cuda_backward(at::Tensor xyz1, at::Tensor xyz2, at::Tensor gradxyz, at::Tensor graddist, at::Tensor idx);
int cuda_auction_find_conflict(at::Tensor xyz1, at::Tensor xyz2, at::Tensor dist, at::Tensor assignment, at::Tensor price,
							   at::Tensor assignment_inv, at::Tensor bid, at::Tensor bid_increments, at::Tensor max_increments,
							   at::Tensor unass_idx, at::Tensor unass_cnt, at::Tensor unass_cnt_sum, at::Tensor cnt_tmp, at::Tensor max_idx, at::Tensor is_conflict, float eps, int iters);
void sia_emd(at::Tensor assignment, at::Tensor assignment_inv, at::Tensor price, at::Tensor xyz1, at::Tensor xyz2, at::Tensor batch_is_conflict, int shutdown_time_milliseconds, long long batch_size, long long n, at::Tensor CURLOOP, at::Tensor number_of_conflict)
{
	int *sia_assignment = (int *)(assignment.data_ptr());
	int *sia_assignment_inv = (int *)(assignment_inv.data_ptr());
	float *pre_p = (float *)(xyz1.data_ptr());
	char *gt_p = (char *)(xyz2.data_ptr());
	float *sia_price = (float *)(price.data_ptr());
	int *batch_is_conflict_p = (int *)(batch_is_conflict.data_ptr());
	int *curloop_ptr = (int *)(CURLOOP.data_ptr());
	int *number_of_conflict_ptr = (int *)(number_of_conflict.data_ptr());
	
	RTree *da[batch_size];
	for (auto i = 0; i < batch_size; i++)
	{

		da[i] = nullptr;
	}

	omp_set_num_threads(batch_size);
    #pragma omp parallel for

	for (auto j = 0; j < batch_size; j++)
	{

		int *assignment_i = &sia_assignment[j * n];
		int *assignment_i_inv = &sia_assignment_inv[j * n];
		char *GT = &gt_p[j * 50];
		float *PRE = &pre_p[j * n * 3];
		float *price_i = &sia_price[j * n];
		int *is_conflict = &batch_is_conflict_p[j * n];
		int *single_curloop = &curloop_ptr[j];
		int *single_conflict = &number_of_conflict_ptr[j];
		int qrysize, objsize;
		CSpatialCCA test;
		test.batch_size = batch_size;
		test.item_size = n;
		qrysize = objsize = n;

		// const char *qrytype, *objtype, *processtype;
		// , qryk, objk;
		// int output;
		// int no_qryupdates, no_objupdates;
		// float randomized_ratio, update_ratio, insert_ratio;

		// ConfigType cr;
		// AddConfigFromFile(cr, "config.ini");
		// AddConfigFromCmdLine(cr, argc, argv);

		// output = getConfigInt("output", cr);
		//	qrytype = getConfigStr("qrytype", cr);
		// objtype = getConfigStr("objtype", cr);
		// qrysize = getConfigInt("qrysize", cr);
		//	objsize = getConfigInt("objsize", cr);
		// qryk = getConfigInt("qryk", cr);
		// objk = getConfigInt("objk", cr);
		// processtype = getConfigStr("process", cr);

		// int dom = getConfigInt("dom", cr);

		// no_qryupdates = getConfigInt("no_qryupdates", cr);
		//	no_objupdates = getConfigInt("no_objupdates", cr);
		// randomized_ratio = getConfigFloat("randomized_ratio", cr);
		// update_ratio = getConfigFloat("update_ratio", cr);
		//	insert_ratio = getConfigFloat("insert_ratio", cr);

		char assignfile[255];
		double start_secs;

		test.readObj(da[j], GT, objsize);
		test.readQry(qrysize, objsize, PRE);
		int ca = 0, cb = 0;

		for (int i = 0; i < qrysize; i++)
		{
			capQ[i] = 1;
		}
		ca = qrysize; // ca += capQ[i]; in the loop
		for (int j = 0; j < objsize; j++)
		{
			capP[j] = 1;
		}
		cb = objsize;
		start_secs = gettime();
		test.initSIA(SIA, qrysize, objsize, capQ, capP);
		test.shutdown_time_milliseconds = shutdown_time_milliseconds;
		auto _ = test.runSIA(is_conflict,single_curloop,0);
		double res = IntToFloat(test.IPA(assignment_i, assignment_i_inv, price_i, qrysize));
		delete da[j];
	}
	
	return;
}

int emd_forward(at::Tensor xyz1, at::Tensor xyz2, at::Tensor dist, at::Tensor assignment, at::Tensor price,
				at::Tensor assignment_inv, at::Tensor bid, at::Tensor bid_increments, at::Tensor max_increments,
				at::Tensor unass_idx, at::Tensor unass_cnt, at::Tensor unass_cnt_sum, at::Tensor cnt_tmp, at::Tensor max_idx, float eps, int iters)
{
	return emd_cuda_forward(xyz1, xyz2, dist, assignment, price, assignment_inv, bid, bid_increments, max_increments, unass_idx, unass_cnt, unass_cnt_sum, cnt_tmp, max_idx, eps, iters);
}

int emd_backward(at::Tensor xyz1, at::Tensor xyz2, at::Tensor gradxyz, at::Tensor graddist, at::Tensor idx)
{

	return emd_cuda_backward(xyz1, xyz2, gradxyz, graddist, idx);
}

int auction_find_conflict(at::Tensor xyz1, at::Tensor xyz2, at::Tensor dist, at::Tensor assignment, at::Tensor price,
						  at::Tensor assignment_inv, at::Tensor bid, at::Tensor bid_increments, at::Tensor max_increments,
						  at::Tensor unass_idx, at::Tensor unass_cnt, at::Tensor unass_cnt_sum, at::Tensor cnt_tmp, at::Tensor max_idx, at::Tensor is_conflict, float eps, int iters)
{
	return cuda_auction_find_conflict(xyz1, xyz2, dist, assignment, price, assignment_inv, bid, bid_increments, max_increments, unass_idx, unass_cnt, unass_cnt_sum, cnt_tmp, max_idx, is_conflict, eps, iters);
}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m)
{
	m.def("find_conflict", &auction_find_conflict, "auction_find_conflict(CUDA)");
	m.def("sia_emd", &sia_emd, "sia");
	m.def("forward", &emd_forward, "emd forward (CUDA)");
	m.def("backward", &emd_backward, "emd backward (CUDA)");
}
