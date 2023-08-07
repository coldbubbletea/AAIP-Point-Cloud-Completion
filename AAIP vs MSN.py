import os 
import numpy as np
import sys
import torch
import time 
sys.path.append('./emd/')
import algorithms
# import copy
# from metrics.metric import l1_cd, l2_cd, emd, f_score
# from Hungarian import *
# from metrics.loss import cd_loss_L1, emd_loss


if __name__ == '__main__':
    lst_txt = './data/compare.list'
    npoints = 1024
    time_limit_microseconds = 20000
    with open(lst_txt, 'r') as f_name:
        name_lst = f_name.readlines()
    
    
    sum_iteration = 0 
    sum_auc_time = 0
    cnt = 0
    log_path = './fortest_%d_%d.txt' % (npoints, time_limit_microseconds)
    log_time_path = './100fortest_time_pre_%d_compare_%d_wo_findconflict.txt' % (npoints, time_limit_microseconds)
    log_aaip_auction_path = './100fortest_time_auc_in_aaip_%d.txt' % npoints
    with open(log_path, 'a') as f_log:
        f_log.write('id ------ GIPA_EMD ------ MSN_EMD ------ PCN_EMD ------ Exact_EMD\n')
    with open(log_time_path, 'a') as f_log:
        f_log.write('TIME\nid ------ GIPA_EMD ------ MSN_EMD ------ PCN_EMD\n')
    total_gipa_time, total_msn_time, total_pcn_time = 0, 0, 0
    for name in name_lst: 
        gt_path = os.path.join('./data/GT/gt_%d' % (npoints), name.strip() + '.txt')
        gt = torch.from_numpy(np.loadtxt(gt_path)).cuda().unsqueeze(0) 
        for i in range(8):
            cnt += 1
            pre_path = os.path.join('./data/PRE/pre_%d' % (npoints), name.strip() + '_%d.txt' % i)
            pre = torch.from_numpy(np.loadtxt(pre_path)).cuda().unsqueeze(0)
            id = tuple([name.strip()])
            print("Object ID:",end=" ")
            print(id)
           # print(pre.shape, gt.shape)
            GIPA = algorithms.SiAuction()
            torch.cuda.synchronize() 
            msn_st = time.time()
            dist = GIPA.Auc(pre, gt, 0.005, 50) # Auc-emd
            emd_msn = torch.sqrt(dist).mean(1).mean().item()
            torch.cuda.synchronize() 
            msn_time = time.time() - msn_st
            print('[EMD/MSN]: %.8f' % emd_msn)
            torch.cuda.synchronize() 
            gipa_st = time.time()
            dist, _, _, auc_iter, auc_time = GIPA.SpAI(pre, gt, id, 0.005, 700, time_limit_microseconds) # gipa-emd
            torch.cuda.synchronize() 
            gipa_time = time.time() - gipa_st
            with open(log_aaip_auction_path, 'a') as f_auc:
                f_auc.write('%s %.6f\n' % (name.strip(), auc_time))
            sum_iteration += auc_iter
            sum_auc_time += auc_time
            emd_gipa = torch.sqrt(dist).mean(1).mean().item()
            
            print('[EMD/AAIP]: %.8f' % emd_gipa)
            # print(emd_msn, emd_gipa)
            # pcn_st = time.time()
            # emd_pcn = emd_loss(pre, gt) # return torch.size([])
            emd_pcn = -1.0
            # pcn_time = time.time() - pcn_st
            pcn_time = -1.0
            # print('[EMD/PCN]: %.8f\nStarting Calculate Exact EMD...' % emd_pcn)
            # matrix, cost_matrix_time = costMatrix(np.array(gt.squeeze(0).cpu()), np.array(pre.squeeze(0).cpu()))
            # print('Finished Cost Matrix...')
            # emd_exact = h_assign(matrix)
            emd_exact = -1.0
            # print('[EMD/Exact]: %.8f' % emd_exact)
            total_gipa_time += gipa_time
            total_msn_time += msn_time
            total_pcn_time += pcn_time
            with open(log_path, 'a') as f_log:
               id_str = name.strip()
               report = '%s %.8f %.8f %.8f %.8f\n' % (id_str, emd_gipa, emd_msn, emd_pcn, emd_exact)
               f_log.write(report)     
            with open(log_time_path, 'a') as f_log_time:
               id_str = name.strip() 
               report = '%s %.8f %.8f %.8f\n' % (id_str, gipa_time, msn_time, pcn_time)
               f_log_time.write(report) 
    with open(log_time_path, 'a') as f_log_time:
        f_log_time.write('Avg. time: %.4f %.4f %.4f\n' % (total_gipa_time / cnt, total_msn_time / cnt, total_pcn_time / cnt))
    print(sum_iteration / 960) 
    
    with open(log_aaip_auction_path, 'a') as f_auc:
        f_auc.write('Avg. iters:[%.2f] time:[%.6f]\n%d cases.\n' % (sum_iteration / cnt, sum_auc_time / cnt, cnt))

