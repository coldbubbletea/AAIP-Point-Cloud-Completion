import time
import numpy as np
import torch
from torch import nn
from torch.autograd import *
import emd
import scipy

class SiAuction():  

    def sia(self, xyz1, xyz2, batch_is_conflict, shutdown_time_milliseconds):
        batchsize, n, _ = xyz1.shape
        _, m, _ = xyz2.shape
        assert(n == m)
        sia_num = torch.zeros(batchsize, device='cpu', dtype=torch.int32)
        conflict_num = torch.zeros(batchsize, device='cpu', dtype=torch.int32)
        price = torch.zeros(batchsize, m, device='cpu').contiguous()
        assignment = torch.zeros(batchsize, m, device='cpu', dtype=torch.int32).contiguous() - 1
        assignment_inv = torch.zeros(batchsize, m, device='cpu', dtype=torch.int32).contiguous() - 1
        emd.sia_emd(assignment, assignment_inv, price, xyz1, xyz2, batch_is_conflict, shutdown_time_milliseconds, batchsize, n, sia_num, conflict_num)
        return assignment, assignment_inv, price, sia_num, conflict_num

    def find_conflict(self, xyz1, xyz2):
        batchsize, n, _ = xyz1.size()
        _, m, _ = xyz2.size()

        assert(n == m)
        assert(xyz1.size()[0] == xyz2.size()[0])
        assert(n % 1024 == 0)
        assert(batchsize <= 512)

        xyz1 = xyz1.contiguous().float().cuda()
        xyz2 = xyz2.contiguous().float().cuda()
        dist = torch.zeros(batchsize, n, device='cuda').contiguous()
        bid = torch.zeros(batchsize, n, device='cuda', dtype=torch.int32).contiguous()
        bid_increments = torch.zeros(batchsize, n, device='cuda').contiguous()
        max_increments = torch.zeros(batchsize, m, device='cuda').contiguous()
        unass_idx = torch.zeros(batchsize * n, device='cuda', dtype=torch.int32).contiguous()
        max_idx = torch.zeros(batchsize * m, device='cuda', dtype=torch.int32).contiguous()
        unass_cnt = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()
        unass_cnt_sum = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()
        cnt_tmp = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()  
        assignment = torch.zeros(batchsize, m, device='cuda', dtype=torch.int32).contiguous() - 1   
        price = torch.zeros(batchsize, m, device='cuda').contiguous()
        assignment_inv = torch.zeros(batchsize, m, device='cuda', dtype=torch.int32).contiguous() - 1
        is_conflict = torch.zeros(batchsize, m, device='cuda', dtype=torch.int32).contiguous()

        emd.find_conflict(xyz1, xyz2, dist, assignment, price, assignment_inv, bid, bid_increments, max_increments, unass_idx, unass_cnt, unass_cnt_sum, cnt_tmp, max_idx, is_conflict, 0.005, 50)
        return is_conflict


class emdFunction(Function):
    @staticmethod
    def forward(ctx, xyz1, xyz2, assignment, assignment_inv, price, eps, iters):

        batchsize, n, _ = xyz1.size()
        _, m, _ = xyz2.size()

        assert(n == m)
        assert(xyz1.size()[0] == xyz2.size()[0])
        assert(n % 1024 == 0)
        assert(batchsize <= 512)

        xyz1 = xyz1.contiguous().float().cuda()
        xyz2 = xyz2.contiguous().float().cuda()
        dist = torch.zeros(batchsize, n, device='cuda').contiguous()
        bid = torch.zeros(batchsize, n, device='cuda', dtype=torch.int32).contiguous()
        bid_increments = torch.zeros(batchsize, n, device='cuda').contiguous()
        max_increments = torch.zeros(batchsize, m, device='cuda').contiguous()
        unass_idx = torch.zeros(batchsize * n, device='cuda', dtype=torch.int32).contiguous()
        max_idx = torch.zeros(batchsize * m, device='cuda', dtype=torch.int32).contiguous()
        unass_cnt = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()
        unass_cnt_sum = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()
        cnt_tmp = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()
        
        emd.forward(xyz1, xyz2, dist, assignment, price, assignment_inv, bid, bid_increments, max_increments, unass_idx, unass_cnt, unass_cnt_sum, cnt_tmp, max_idx, eps, iters)

        ctx.save_for_backward(xyz1, xyz2, assignment)
        return dist, assignment

    @staticmethod
    def backward(ctx, graddist, gradidx):
        xyz1, xyz2, assignment = ctx.saved_tensors
        graddist = graddist.contiguous()

        gradxyz1 = torch.zeros(xyz1.size(), device='cuda').contiguous()
        gradxyz2 = torch.zeros(xyz2.size(), device='cuda').contiguous()

        emd.backward(xyz1, xyz2, gradxyz1, graddist, assignment)
        return gradxyz1, gradxyz2, None, None, None, None, None



class emdModule(nn.Module):
    def __init__(self):
        super(emdModule, self).__init__()

    def forward(self, input1, input2, ass, ass_inv, price, eps, iters):
        return emdFunction.apply(input1, input2, ass, ass_inv, price, eps, iters)

def test_emd():
    xyz1 = torch.rand(32, 8192, 3, device='cuda', requires_grad=True)
    xyz2 = torch.rand(32, 8192, 3, device='cuda', requires_grad=False)
    # print(xyz1.shape, xyz2.shape)
    # print(xyz1[0, :10, :], xyz2[0, :10, :])
    # assignment,assignment_inv,price = sia(x1.cpu(), x2.cpu(), 10)
    # xyz1 = torch.from_numpy(np.loadtxt('/home/steven/Workspace/PointCloudCompletion_SIA-auction/rank_output/rank03/network_0/predicted_0.txt')).unsqueeze(0)
    # xyz2 = torch.from_numpy(np.loadtxt('/home/steven/Workspace/PointCloudCompletion_SIA-auction/rank_output/rank03/gt/gt_0.txt')).unsqueeze(0)
    # print(xyz1.shape, xyz2.shape)
    # print(xyz1[0, :10, :], xyz2[0, :10, :])
    SiA = SiAuction(xyz1, xyz2)
    batchsize, n, _ = xyz1.size()
    _, m, _ = xyz2.size()
    assert(n == m)
    assert(xyz1.size()[0] == xyz2.size()[0])
    assert(n % 1024 == 0)
    assert(batchsize <= 512)
    auction_emd = emdModule()
    xyz1 = xyz1.contiguous().float().cuda()
    xyz2 = xyz2.contiguous().float().cuda()
    dist = torch.zeros(batchsize, n, device='cuda').contiguous()
    bid = torch.zeros(batchsize, n, device='cuda', dtype=torch.int32).contiguous()
    bid_increments = torch.zeros(batchsize, n, device='cuda').contiguous()
    max_increments = torch.zeros(batchsize, m, device='cuda').contiguous()
    assignment = torch.zeros(batchsize, m, device='cuda', dtype=torch.int32).contiguous() - 1
    # is_conflict = torch.zeros(batchsize, m, device='cuda', dtype=torch.int32).contiguous()
    price = torch.zeros(batchsize, m, device='cuda').contiguous()
    assignment_inv = torch.zeros(batchsize, m, device='cuda', dtype=torch.int32).contiguous() - 1
    unass_idx = torch.zeros(batchsize * n, device='cuda', dtype=torch.int32).contiguous()
    max_idx = torch.zeros(batchsize * m, device='cuda', dtype=torch.int32).contiguous()
    unass_cnt = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()
    unass_cnt_sum = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()
    cnt_tmp = torch.zeros(512, dtype=torch.int32, device='cuda').contiguous()    
    # torch.set_printoptions(profile="full")
   
    is_conflict = SiA.find_conflict(xyz1, xyz2)
    print('Conflict_List:{}'.format(is_conflict))
    print(xyz1.size())
    assignment, assignment_inv, price, sia_num, conflict_num = SiA.sia(xyz1.cpu(), xyz2.cpu(), is_conflict.cpu(), 1000)
    print('SIA_Assign:{}'.format(assignment))
    # emdModule(xyz1, xyz2, assignment, assignment_inv, price, 0.005, 50)
    # print(assignment)
    # print(is_conflict)

    # assignment = torch.tensor(assignment, device='cuda')
    # assignment_inv = torch.tensor(assignment_inv, device='cuda')
    # price = torch.tensor(price, device='cuda')
    # assignment = torch.zeros(32, 8192, device='cuda', dtype=torch.int32).contiguous() - 1
    # assignment_inv = torch.zeros(32, 8192, device='cuda', dtype=torch.int32).contiguous() - 1
    # price = torch.zeros(32, 8192, device='cuda', dtype=torch.float32).contiguous()
    dist, _ = auction_emd(xyz1, xyz2, assignment, assignment_inv, price ,0.005, 50) 
    emd1 = (torch.sqrt(dist).mean(1)).cpu().detach().numpy()
    print(emd1)

    # Pure EMD
    # assign and inv -1
    # price 0
    assignment = torch.zeros(32, 8192, device='cuda', dtype=torch.int32).contiguous() - 1
    assignment_inv = torch.zeros(32, 8192, device='cuda', dtype=torch.int32).contiguous() - 1
    price = torch.zeros(32, 8192, device='cuda', dtype=torch.float32).contiguous()
    dist, _ = auction_emd(xyz1,xyz2,assignment, assignment_inv, price ,0.005, 50) 
    emd_auc_approximate = (torch.sqrt(dist).mean(1)).cpu().detach().numpy()
    print(emd_auc_approximate)

    assignment = torch.zeros(32, 8192, device='cuda', dtype=torch.int32).contiguous() - 1
    assignment_inv = torch.zeros(32, 8192, device='cuda', dtype=torch.int32).contiguous() - 1
    price = torch.zeros(32, 8192, device='cuda', dtype=torch.float32).contiguous()
    dist, _ = auction_emd(xyz1, xyz2, assignment, assignment_inv, price, 0.005, 20000) 
    emd_exact = (torch.sqrt(dist).mean(1)).cpu().detach().numpy()
    print(emd_exact)

    rho_1 = scipy.stats.spearmanr(emd_exact, emd_auc_approximate)[0]
    rho_2 = scipy.stats.spearmanr(emd_exact, emd1)[0]
    print('rho_1:{}; rho_2:{}'.format(rho_1, rho_2))

if __name__ == '__main__':
    test_emd()
