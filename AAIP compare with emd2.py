import os 
import numpy as np
import sys
import torch
import time 
sys.path.append('./emd/')
import algorithms


if __name__ == '__main__':
    print("============================================================================")
    print("In this program,\nwe compare the time cost and MSE between AAIP(our method) and emd2.")
    print("Two input point clouds are predict point cloud we grabbed in training and groud point cloud.")
    print("============================================================================")
    lst_txt = './data/compare.list'
    npoints = 8192
    sia_time_control = 200000 #Time to control SIA
    with open(lst_txt, 'r') as f_name:
        name_lst = f_name.readlines()
    sum_iteration = 0 
    sum_auc_time = 0
    cnt = 0
    emd_value_comparison_path = 'emd_value_comparison.txt'
    time_cost_comparison_path = 'time_cost_comparison.txt'
    with open(emd_value_comparison_path, 'a') as f_log:
        f_log.write('emd value\nid ------ AAIP_emd ------ emd2\n')
    with open(time_cost_comparison_path, 'a') as f_log:
        f_log.write('TIME(s)\n id ------ AAIP_emd_time ------ emd2_time\n')
    total_emd_approximation_time, total_msn_time, total_pcn_time = 0, 0, 0
    for name in name_lst: 
        gt_path = os.path.join('./data/GT/gt_%d' % (npoints), name.strip() + '.txt')
        gt = torch.from_numpy(np.loadtxt(gt_path)).cuda().unsqueeze(0) 
        cnt=0
        for i in range(8):
            cnt += 1
            pre_path = os.path.join('./data/PRE/pre_%d' % (npoints), name.strip() + '_%d.txt' % i)
            pre = torch.from_numpy(np.loadtxt(pre_path)).cuda().unsqueeze(0)
            id = tuple([name.strip()])
            emd_approximation = algorithms.SiAuction()
            torch.cuda.synchronize() 
            msn_st = time.time()
            dist = emd_approximation.Auc(pre, gt, 0.005) # Auc-emd
            emd_msn = torch.sqrt(dist).mean(1).mean().item()
            torch.cuda.synchronize() 
            msn_time = time.time() - msn_st
            torch.cuda.synchronize() 
            emd_approximation_st = time.time()
            dist, _, _, auc_iter, auc_time = emd_approximation.SpAI(pre, gt, id, 0.005, 50, sia_time_control) # emd_approximation-emd
            torch.cuda.synchronize() 
            emd_approximation_time = time.time() - emd_approximation_st
            sum_iteration += auc_iter
            sum_auc_time += auc_time
            emd_emd_approximation = torch.sqrt(dist).mean(1).mean().item()
            emd_pcn = -1.0
            pcn_time = -1.0
            emd_exact = -1.0
            total_emd_approximation_time += emd_approximation_time
            total_msn_time += msn_time
            total_pcn_time += pcn_time
            with open(emd_value_comparison_path, 'a') as f_log:
               id_str = name.strip()
               report = '%s %.8f %.8f\n' % (id_str, emd_emd_approximation, emd_msn)
               f_log.write(report)     
            with open(time_cost_comparison_path, 'a') as f_log_time:
               id_str = name.strip() 
               report = '%s %.8f %.8f \n' % (id_str, emd_approximation_time, msn_time)
               f_log_time.write(report) 
            print('object_id: %s %d/8' %(name.strip(),cnt))
    with open(time_cost_comparison_path, 'a') as f_log_time:
        f_log_time.write('Avg. time: %.4f %.4f \n' % (total_emd_approximation_time / 80, total_msn_time / 80))
   