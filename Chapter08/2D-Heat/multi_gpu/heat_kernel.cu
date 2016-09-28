/*
******************************************************
This file is the OpenACC multi-GPU version of 2D Heat Equation 
using OpenMP+OpenACC hybrid model. This implementation is based
on the CPU version from
http://www.many-core.group.cam.ac.uk/archive/CUDAcourse09/

Permission to use, copy, distribute and modify this software for any 
purpose with or without fee is hereby granted. This software is        
provided "as is" without express or implied warranty. 

Send comments or suggestions for this OpenACC version to
            rxu6@uh.edu, schandra@udel.edu

Authors: Rengan Xu, Sunita Chandrasekaran

May 26th, 2016
******************************************************
*/

extern "C" __global__ void step_kernel(int ni, 
                     int nj,
                     double tfac, 
                     double *temp_in,
                     double *temp_out) 
{
    int i, j, i00, im10, ip10, i0m1, i0p1;
    double d2tdx2, d2tdy2;

    j = blockIdx.y + 1;
    while(j < nj-1)
    {
        i = threadIdx.x + blockIdx.x*blockDim.x + 1;
        while(i < ni-1)
        {
            i00 = i + ni*j;
            im10 = i-1 + ni*j;
            ip10 = i+1 + ni*j;
            i0m1 = i + ni*(j-1);
            i0p1 = i + ni*(j+1);

            d2tdx2 = temp_in[im10] - 2*temp_in[i00] + temp_in[ip10];
            d2tdy2 = temp_in[i0m1] - 2*temp_in[i00] + temp_in[i0p1];
            
            temp_out[i00] = temp_in[i00] + tfac*(d2tdx2 + d2tdy2);
            i += blockDim.x*gridDim.x;
        }
        j += gridDim.y;
    }
}
