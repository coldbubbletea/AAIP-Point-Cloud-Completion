import time
import numpy as np
import torch
from torch import nn
from torch.autograd import Function
import emd
def sia(xyz1, xyz2, loop):
    batchsize, n, _ = xyz1.size()
    assignment = torch.zeros(batchsize, n // 2, device='cpu', dtype=torch.int32).contiguous() 
    emd_value = torch.zeros(batchsize,device='cpu', dtype=torch.float32).contiguous() 
    unassignment = emd.SIA_emd(assignment,emd_value,xyz1,xyz2,loop,n)
    # print(assignment.shape)
    # print(unassignment)
    # print(emd_value)
    return assignment,unassignment,emd_value

def assign_points(xyz2, assignment):
    assignment = torch.tensor(assignment, device='cuda')
    bs, n = assignment.shape
    x2 = torch.zeros(assignment.shape[0], assignment.shape[1], 3, device='cuda')
    for i in range(bs):
        x2[i, :, :] = xyz2[i, assignment[i, :].long(), :]
    return x2

def cal_dist(sia_xyz1, sia_xyz2):
    assert sia_xyz1.shape == sia_xyz2.shape
    # batchsize, n, _ = sia_xyz1.shape
    d = torch.sqrt(((sia_xyz1 - sia_xyz2) ** 2).sum(-1)).mean(-1)
    return d

         