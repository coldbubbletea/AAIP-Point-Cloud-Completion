/***********************************************************************************************
 ***                                AAIP       A L G O R I T H M                             ***
 ***********************************************************************************************
 *                                                                                             *
 *                   Paper Name :                                                              *
 *                                                                                             *
 *                   Authors :  Shunran ZHANG, Xiubo ZHANG                                     *
 *                                                                                             *
 *                   File Name : emd_cuda.cu                                                   *
 *                                                                                             *
 *                   Programmer : Xiubo ZHANG                                                  *
 *                                                                                             *
 *                   Start Date : Mar 3, 2023                                                  *
 *                                                                                             *
 *                   Last Update : July 27, 2023                                              *
 *                                                                                             *
 *---------------------------------------------------------------------------------------------*/
#include <stdio.h>
#include <ATen/ATen.h>
#include <cuda.h>
#include <iostream>
#include <cuda_runtime.h>

__device__ __forceinline__ float atomicMax(float *address, float val)
{
    int ret = __float_as_int(*address);
    while(val > __int_as_float(ret))
    {
        int old = ret;
        if((ret = atomicCAS((int *)address, old, __float_as_int(val))) == old)
            break;
    }
    return __int_as_float(ret);
}


__global__ void clear_here(int b, int * cnt_tmp, int * unass_cnt) {
	for (int i = threadIdx.x; i < b; i += blockDim.x) {
		cnt_tmp[i] = 0;
		unass_cnt[i] = 0;
	}
}

__global__ void calc_unass_cnt_here(int b, int n, int * assignment, int * unass_cnt) { 
	// count the number of unassigned points in each batch
	const int BLOCK_SIZE = 1024; 
	__shared__ int scan_array[BLOCK_SIZE];
	for (int i = blockIdx.x; i < b; i += gridDim.x) {
		scan_array[threadIdx.x] = assignment[i * n + blockIdx.y * BLOCK_SIZE + threadIdx.x] == -1 ? 1 : 0;
		__syncthreads();
		
		int stride = 1;
		while(stride <= BLOCK_SIZE / 2) {
			int index = (threadIdx.x + 1) * stride * 2 - 1; 
			if(index < BLOCK_SIZE)
				scan_array[index] += scan_array[index - stride]; 
			stride = stride * 2;
			__syncthreads(); 
		}
		__syncthreads();
		
		if (threadIdx.x == BLOCK_SIZE - 1) {
			atomicAdd(&unass_cnt[i], scan_array[threadIdx.x]);
		}
		__syncthreads();
	}
}

__global__ void calc_unass_cnt_sum_here(int b, int * unass_cnt, int * unass_cnt_sum) {
	// count the cumulative sum over over unass_cnt
	const int BLOCK_SIZE = 512; // batch_size <= 512
	__shared__ int scan_array[BLOCK_SIZE];
	scan_array[threadIdx.x] = unass_cnt[threadIdx.x];
	__syncthreads();
	
	int stride = 1;
	while(stride <= BLOCK_SIZE / 2) {
		int index = (threadIdx.x + 1) * stride * 2 - 1; 
		if(index < BLOCK_SIZE)
			scan_array[index] += scan_array[index - stride]; 
		stride = stride * 2;
		__syncthreads(); 
	}
	__syncthreads();
	stride = BLOCK_SIZE / 4; 
	while(stride > 0) {
		int index = (threadIdx.x + 1) * stride * 2 - 1; 
		if((index + stride) < BLOCK_SIZE)
			scan_array[index + stride] += scan_array[index];
		stride = stride / 2;
		__syncthreads(); 
	}
	__syncthreads(); 
	
	//printf("%d\n", unass_cnt_sum[b - 1]);
	unass_cnt_sum[threadIdx.x] = scan_array[threadIdx.x];
}

__global__ void calc_unass_idx_here(int b, int n, int * assignment, int * unass_idx, int * unass_cnt, int * unass_cnt_sum, int * cnt_tmp) {
	// list all the unassigned points
	for (int i = blockIdx.x; i < b; i += gridDim.x) {
		if (assignment[i * n + blockIdx.y * 1024 + threadIdx.x] == -1) {
			int idx = atomicAdd(&cnt_tmp[i], 1);
			unass_idx[unass_cnt_sum[i] - unass_cnt[i] + idx] = blockIdx.y * 1024 + threadIdx.x;
		} 
	}
}

__global__ void Bid_here(int b, int n, const float * xyz1, const float * xyz2, float eps, int * assignment, int * assignment_inv, float * price, 
					int * bid, float * bid_increments, float * max_increments, int * unass_cnt, int * unass_cnt_sum, int * unass_idx) {
	const int batch = 2048, block_size = 1024, block_cnt = n / 1024;
	__shared__ float xyz2_buf[batch * 3];
	__shared__ float price_buf[batch];
	__shared__ float best_buf[block_size];
	__shared__ float better_buf[block_size];
	__shared__ int best_i_buf[block_size];
	for (int i = blockIdx.x; i < b; i += gridDim.x) {
		int _unass_cnt = unass_cnt[i];
		if (_unass_cnt == 0)
			continue;
		int _unass_cnt_sum = unass_cnt_sum[i];
		int unass_per_block = (_unass_cnt + block_cnt - 1) / block_cnt;
		int thread_per_unass = block_size / unass_per_block;
		int unass_this_block = max(min(_unass_cnt - (int) blockIdx.y * unass_per_block, unass_per_block), 0);
			
		float x1, y1, z1, best = -1e9, better = -1e9;
		int best_i = -1, _unass_id = -1, thread_in_unass;

		if (threadIdx.x < thread_per_unass * unass_this_block) {
			_unass_id = unass_per_block * blockIdx.y + threadIdx.x / thread_per_unass + _unass_cnt_sum - _unass_cnt;
			_unass_id = unass_idx[_unass_id];
			thread_in_unass = threadIdx.x % thread_per_unass;

			x1 = xyz1[(i * n + _unass_id) * 3 + 0];
			y1 = xyz1[(i * n + _unass_id) * 3 + 1];
			z1 = xyz1[(i * n + _unass_id) * 3 + 2];
		}

		for (int k2 = 0; k2 < n; k2 += batch) {
			int end_k = min(n, k2 + batch) - k2;
			for (int j = threadIdx.x; j < end_k * 3; j += blockDim.x) {
				xyz2_buf[j] = xyz2[(i * n + k2) * 3 + j];
			}
			for (int j = threadIdx.x; j < end_k; j += blockDim.x) {
				price_buf[j] = price[i * n + k2 + j];
			}
			__syncthreads();

			if (_unass_id != -1) {
				int delta = (end_k + thread_per_unass - 1) / thread_per_unass;
				int l = thread_in_unass * delta;
				int r = min((thread_in_unass + 1) * delta, end_k);
				for (int k = l; k < r; k++) 
				//if (!last || assignment_inv[i * n + k + k2] == -1)
				{
					float x2 = xyz2_buf[k * 3 + 0] - x1;
					float y2 = xyz2_buf[k * 3 + 1] - y1;
					float z2 = xyz2_buf[k * 3 + 2] - z1;
					// the coordinates of points should be normalized to [0, 1]
					float d = 3.0 - sqrtf(x2 * x2 + y2 * y2 + z2 * z2) - price_buf[k];
					if (d > best) {
						better = best;
						best = d;
						best_i = k + k2;
					}
					else if (d > better) {
						better = d;
					}
				}
			}
			__syncthreads();
		}

		best_buf[threadIdx.x] = best;
		better_buf[threadIdx.x] = better;
		best_i_buf[threadIdx.x] = best_i;
		__syncthreads();
		
		if (_unass_id != -1 && thread_in_unass == 0) {
			for (int j = threadIdx.x + 1; j < threadIdx.x + thread_per_unass; j++) {
				if (best_buf[j] > best) {
					better = max(best, better_buf[j]);
					best = best_buf[j];
					best_i = best_i_buf[j];
				}
				else better = max(better, best_buf[j]);
			}
			bid[i * n + _unass_id] = best_i;
			bid_increments[i * n + _unass_id] = best - better + eps; 
			atomicMax(&max_increments[i * n + best_i], best - better + eps);
		}
	}
}





int cuda_auction_find_conflict(at::Tensor xyz1, at::Tensor xyz2, at::Tensor dist, at::Tensor assignment, at::Tensor price, 
	                 at::Tensor assignment_inv, at::Tensor bid, at::Tensor bid_increments, at::Tensor max_increments,
	                 at::Tensor unass_idx, at::Tensor unass_cnt, at::Tensor unass_cnt_sum, at::Tensor cnt_tmp, at::Tensor max_idx,at::Tensor is_conflict, float eps, int iters) {

	const auto batch_size = xyz1.size(0);
	const auto n = xyz1.size(1); //num_points point cloud A
	const auto m = xyz2.size(1); //num_points point cloud B
	
	if (n != m) {
		printf("Input Error! The two point clouds should have the same size.\n");
		return -1;
	}

	if (batch_size > 512) {
		printf("Input Error! The batch size should be less than 512.\n");
		return -1;
	}

	if (n % 1024 != 0) {
		printf("Input Error! The size of the point clouds should be a multiple of 1024.\n");
		return -1;
	}

	//cudaEvent_t start,stop;
	//cudaEventCreate(&start);
	//cudaEventCreate(&stop);
	//cudaEventRecord(start);
	//int iters = 50;
	// SIA_PRE_ASSIGN<<<dim3(batch_size, n / 1024, 1), 1024>>>(batch_size, n, xyz1.data<float>(), xyz2.data<float>(), eps, assignment.data<int>(), assignment_inv.data<int>(), 
	// 		                          price.data<float>(), bid.data<int>(), bid_increments.data<float>(), max_increments.data<float>(),
	// 		                          unass_cnt.data<int>(), unass_cnt_sum.data<int>(), unass_idx.data<int>());
        
		clear_here<<<1, batch_size>>>(batch_size, cnt_tmp.data<int>(), unass_cnt.data<int>());
		calc_unass_cnt_here<<<dim3(batch_size, n / 1024, 1), 1024>>>(batch_size, n, assignment.data<int>(), unass_cnt.data<int>());
		calc_unass_cnt_sum_here<<<1, batch_size>>>(batch_size, unass_cnt.data<int>(), unass_cnt_sum.data<int>());
		calc_unass_idx_here<<<dim3(batch_size, n / 1024, 1), 1024>>>(batch_size, n, assignment.data<int>(), unass_idx.data<int>(), unass_cnt.data<int>(), 
											 unass_cnt_sum.data<int>(), cnt_tmp.data<int>());
		Bid_here<<<dim3(batch_size, n / 1024, 1), 1024>>>(batch_size, n, xyz1.data<float>(), xyz2.data<float>(), eps, assignment.data<int>(), assignment_inv.data<int>(), 
			                          price.data<float>(), bid.data<int>(), bid_increments.data<float>(), max_increments.data<float>(),
			                          unass_cnt.data<int>(), unass_cnt_sum.data<int>(), unass_idx.data<int>());
		
		
	

    
	cudaError_t err = cudaGetLastError();
	  if (err != cudaSuccess) {
	    printf("error in nnd Output: %s\n", cudaGetErrorString(err));
	    return 0;
	  }
	  return 1;
}





