from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CUDAExtension

setup(
    name='emd',
    ext_modules=[
        
        CUDAExtension('emd', [
            'emd.cpp',    # Same function as main.cpp
            'SpatialCCA.cpp',
            'OptConstraint.cpp',
            'OptAssign.cpp',
            'nheap.cpp',
            'OptUpdate.cpp',
            'rtree/rtree.cc',
            'rtree/hilbert.cc',
            'rtree/blk_file.cc',
            'rtree/functions.cc',
            'rtree/linlist.cc',           
            'emd_cuda.cu',
            'auction_find_conflict.cu'
           

        ], extra_compile_args={'cxx':['-g','-fopenmp'],'nvcc':['-O3']} ),
    ],
    cmdclass={
        'build_ext': BuildExtension
    })





