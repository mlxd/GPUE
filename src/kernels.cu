/*** kernels.cu - GPUE: Split Operator based GPU solver for Nonlinear 
Schrodinger Equation, Copyright (C) 2011-2015, Lee J. O'Riordan 
<loriordan@gmail.com>, Tadhg Morgan, Neil Crowley. 
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

1. Redistributions of source code must retain the above copyright 
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its 
contributors may be used to endorse or promote products derived from 
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "../include/constants.h"
#include <stdio.h>


//Evaluted in MATLAB: N*4*HBAR*HBAR*PI*(4.67e-9/mass)*sqrt(mass*(omegaZ)/(2*PI*HBAR))
__constant__ double gDenConst = 6.6741e-40;

inline __device__ unsigned int getGid3d3d(){
    return blockDim.x * ( ( blockDim.y * ( ( blockIdx.z * blockDim.z + threadIdx.z ) + blockIdx.y ) + threadIdx.y ) + blockIdx.x ) + threadIdx.x;
}

// function to perform a transposition (2d) or permutation (3d)
// Note: The 3 ints represent the final placement of that data direction
//       after transposition
inline __device__ unsigned int permuteGid(int d1, int d2, int d3){

    // I cannot seem to think of any way to write this in a general case...

    unsigned int x, y, z;

    // If the three axes are in the original directions.
    if (d1 == 0 && d2 == 1 && d3 == 2){
        return getGid3d3d();
    } 

    else if (d1 == 1 && d2 == 2 && d3 == 0){
        x = blockIdx.x * blockDim.x + threadIdx.x;
        z = blockDim.z * (x + blockIdx.z) + threadIdx.z;
        y = blockDim.y * (z + blockIdx.y) + threadIdx.y;
        return y;
    }

    else if (d1 == 2 && d2 == 0 && d3 == 1){
        y = blockIdx.y * blockDim.y + threadIdx.y;
        x = blockDim.x * (y + blockIdx.x) + threadIdx.x;
        z = blockDim.z * (x + blockIdx.z) + threadIdx.z;
        return z;
    }

    else if (d1 == 0 && d2 == 2 && d3 == 1){
        y = blockIdx.y * blockDim.y + threadIdx.y;
        z = blockDim.z * (y + blockIdx.z) + threadIdx.z;
        x = blockDim.x * (z + blockIdx.x) + threadIdx.x;
        return x;
    }

    else if (d1 == 1 && d2 == 0 && d3 == 2){
        //z = blockIdx.z * blockDim.z + threadIdx.z;
        //y = blockDim.y * (z + blockIdx.x) + threadIdx.y;
        //x = blockDim.x * (y + blockIdx.y) + threadIdx.x;
        //x = blockDim.x * (z + blockIdx.x) + threadIdx.x;
        //y = blockDim.y * (x + blockIdx.y) + threadIdx.y;
        //return y;
        //return x;
        z = blockIdx.z*blockDim.z + threadIdx.z;
        x = blockIdx.x*blockDim.x + threadIdx.x;
        y = blockIdx.y*blockDim.y + threadIdx.y;
        return x + blockDim.x*y + blockDim.y*blockDim.x*z;
    }

    else if (d1 == 2 && d2 == 1 && d3 == 0){
        x = blockIdx.x * blockDim.x + threadIdx.x;
        y = blockDim.y * (x + blockIdx.y) + threadIdx.y;
        z = blockDim.z * (y + blockIdx.z) + threadIdx.z;
        return z;
    }
    else{
        return 0;
    }

    
}

__device__ unsigned int getBid3d3d(){
    return blockIdx.x + gridDim.x*(blockIdx.y + gridDim.y * blockIdx.z);
}

__device__ unsigned int getTid3d3d(){
    return blockDim.x * ( blockDim.y * ( blockDim.z + ( threadIdx.z * blockDim.y ) )  + threadIdx.y )  + threadIdx.x;
}

__device__ double2 conjugate(double2 in){
    double2 result = in;
    result.y = -result.y;
    return result;
}

__device__ double2 realCompMult(double scalar, double2 comp){
    double2 result;
    result.x = scalar * comp.x;
    result.y = scalar * comp.y;
    return result;
}

__device__ double complexMagnitude(double2 in){
    return sqrt(in.x*in.x + in.y*in.y);
}

__host__ __device__ double complexMagnitudeSquared(double2 in){
    return in.x*in.x + in.y*in.y;
}

__host__ __device__ double2 complexMultiply(double2 in1, double2 in2){
    double2 result;
    result.x = (in1.x*in2.x - in1.y*in2.y);
    result.y = (in1.x*in2.y + in1.y*in2.x);
    return result;
}

/*
* Used to perform conj(in1)*in2; == < in1 | in2 >
*/
inline __device__ double2 braKetMult(double2 in1, double2 in2){
    return complexMultiply(conjugate(in1),in2);
}

/**
 * Performs complex multiplication of in1 and in2, giving result as out. 
 */
__global__ void cMult(double2* in1, double2* in2, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    double2 tin1 = in1[gid];
    double2 tin2 = in2[gid];
    result.x = (tin1.x*tin2.x - tin1.y*tin2.y);
    result.y = (tin1.x*tin2.y + tin1.y*tin2.x);
    out[gid] = result;
}

__global__ void cMultPhi(double2* in1, double* in2, double2* out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = cos(in2[gid])*in1[gid].x - in1[gid].y*sin(in2[gid]);
    result.y = in1[gid].x*sin(in2[gid]) + in1[gid].y*cos(in2[gid]);
    out[gid] = result;
}

/**
 * Performs multiplication of double* with double2*
 */
__global__ void vecMult(double2 *in, double *factor, double2 *out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = (in[gid].x * factor[gid]);
    result.y = (in[gid].y * factor[gid]);
    out[gid] = result;
}


/**
 * Performs the non-linear evolution term of Gross--Pitaevskii equation.
 */
__global__ void cMultDensity(double2* in1, double2* in2, double2* out, double dt, double mass,double omegaZ, int gstate, int N){
    double2 result;
    double gDensity;
    int gid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
    double2 tin1 = in1[gid];
    double2 tin2 = in2[gid];
    gDensity = gDenConst*complexMagnitudeSquared(in2[gid])*(dt/HBAR); // scaling of interaction strength doesn't work now

    if(gstate == 0){
        double tmp = in1[gid].x*exp(-gDensity);
        result.x = (tmp)*tin2.x - (tin1.y)*tin2.y;
        result.y = (tmp)*tin2.y + (tin1.y)*tin2.x;
    }
    else{
        double2 tmp;
        tmp.x = tin1.x*cos(-gDensity) - tin1.y*sin(-gDensity);
        tmp.y = tin1.y*cos(-gDensity) + tin1.x*sin(-gDensity);
        
        result.x = (tmp.x)*tin2.x - (tmp.y)*tin2.y;
        result.y = (tmp.x)*tin2.y + (tmp.y)*tin2.x;
    }
    out[gid] = result;
}

/**
 * Divides both components of vector type "in", by the value "factor".
 * Results given with "out".
 */
__global__ void scalarDiv(double2* in, double factor, double2* out){
    double2 result;
    //extern __shared__ double2 tmp_in[];
    unsigned int gid = getGid3d3d();
    result.x = (in[gid].x / factor);
    result.y = (in[gid].y / factor);
    out[gid] = result;
}

/**
 * Multiplies both components of vector type "in", by the value "factor".
 * Results given with "out". 
 */
__global__ void scalarMult(double2* in, double factor, double2* out){
    double2 result;
    //extern __shared__ double2 tmp_in[];
    unsigned int gid = getGid3d3d();
    result.x = (in[gid].x * factor);
    result.y = (in[gid].y * factor);
    out[gid] = result;
}

/**
 * As above, but normalises for wfc
 */
__global__ void scalarDiv_wfcNorm(double2* in, double dr, double2* pSum, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    double norm = sqrt((pSum[0].x + pSum[0].y)*dr);
    result.x = (in[gid].x/norm);
    result.y = (in[gid].y/norm);
    out[gid] = result;
}

/**
 * Raises in to the power of param
 */
__global__ void scalarPow(double2* in, double param, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    result.x = pow(result.x, param);
    result.y = pow(result.y, param);
    out[gid] = result;
}

/**
 * Finds conjugate for double2*
 */
__global__ void vecConjugate(double2 *in, double2 *out){
    double2 result;
    unsigned int gid = getGid3d3d(); 
    result.y = -in[gid].y;
    out[gid] = result;
}

__global__ void angularOp(double omega, double dt, double2* wfc, double* xpyypx, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    double op;
    op = exp( -omega*xpyypx[gid]*dt);
    result.x=wfc[gid].x*op;
    result.y=wfc[gid].y*op;
    out[gid]=result;
}

/**
 * Routine for parallel summation. Can be looped over from host.
 */
__global__ void multipass(double2* input, double2* output, int pass){
    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x;
    // printf("bid0=%d\n",bid);

    unsigned int gid = getGid3d3d();
    extern __shared__ double2 sdata[];
    sdata[tid] = input[gid];
    if(pass == 0){
        sdata[tid].x *= sdata[tid].x;
        sdata[tid].y *= sdata[tid].y;    
    }
    __syncthreads();
    for(int i = blockDim.x>>1; i > 0; i>>=1){
        if(tid < blockDim.x>>1){
            sdata[tid].x += sdata[tid + i].x;
            sdata[tid].y += sdata[tid + i].y;
        }
        __syncthreads();
    }
    if(tid==0){
        output[bid] = sdata[0];
    }
}


/*
* Calculates all of the energy of the current state. sqrt_omegaz_mass = sqrt(omegaZ/mass), part of the nonlin interaction term
*/
__global__ void energyCalc(double2 *wfc, double2 *op, double dt, double2 *energy, int gnd_state, int op_space, double sqrt_omegaz_mass){
    unsigned int gid = getGid3d3d();
    double hbar_dt = HBAR/dt;
    double g_local = 0.0;
    double2 result;
    double opLocal;
    if(op_space)
        g_local = gDenConst*sqrt_omegaz_mass*complexMagnitudeSquared(wfc[gid]);
    if(!gnd_state){
        opLocal = -log(op[gid].x + g_local)*hbar_dt;
    }
    else{
        opLocal = cos(op[gid].x + g_local)*hbar_dt;
    }
    result = braKetMult(wfc[gid], realCompMult(opLocal,wfc[gid]));
    //printf("oplocal=%e    Resx=%e    Resy=%e\n",opLocal,result.x,result.y);
    energy[gid].x += result.x;
    energy[gid].y += result.y;
}

// Kernel for 2d transpose, note global for now...
__global__ void transpose2d(double *indata, double *outdata){
    unsigned int gid1 = getGid3d3d();

    // Note that this should always be 102 for 2d transpose if we ignore z
    unsigned int gid2 = permuteGid(1,0,2);

    outdata[gid2] = indata[gid1];
}

// Kernel for 2d transpose, note global for now...
__global__ void naivetranspose2d(int xDim, int yDim, 
                            const double *indata, double *outdata){
    const double* tmp = indata;
    int index = 0;
    int index2 = 0;
    for (int i = 0; i < xDim; i++){
        for (int j = 0; j < yDim; j++){
            index = j + i*yDim;
            index2 = i + j*xDim;
            outdata[index] = tmp[index2];
        }
    }
}

__global__ void transpose2d2(const double2 *indata, double2 *outdata){
    unsigned int gid1 = getGid3d3d();

    // Note that this should always be 102 for 2d transpose if we ignore z
    unsigned int gid2 = permuteGid(0,1,2);

    outdata[gid1] = indata[gid2];
}

__global__ void naivetranspose2d2(int xDim, int yDim, 
                             const double2 *indata, double2 *outdata){
    const double2 *tmp = indata;
    int index = 0;
    int index2 = 0;
    for (int i = 0; i < xDim; i++){
        for (int j = 0; j < yDim; j++){
            index = j + i*yDim;
            index2 = i + j*xDim;
            outdata[index] = tmp[index2];
        }
    }
}

/*
 * Calculates the trapezoidalm integration of a provided double* vector
 * Note: -This is a variation on figure 9.7 on page 209 of Programming Massively
 *        Parallel Processors, Second Edition.
 *       -We do not complete the entire prefix sum because we only need the 
 *        final sum in the end. 
 *       -Further modifications have been made for the trapezoidal rule...
 */
/*
__global__ void trapz(double *array, const int dimension, double dx, 
                      double out){

    // src material uses global variable SECTION_SIZE. dunno...
    double *element = array;

    int i = blockIdx.x*blockDim.x + threadIdx.x;

    // src material uses InputSize here. Not sure what it all means
    if (i < dimension){
        element[threadIdx.x] = array[i];
    }

    for (unsigned int stride = 1; stride < blockIdx.x; stride *= 2){
        __syncthreads();
        int index = (threadIdx.x+1) * 2 * stride - 1;
        if (index < blockDim.x){
            element[index + stride] += element[index]; 
        }
    }

    __syncthreads();

    out = element[dimension];
}
*/

//##############################################################################
//##############################################################################

/**
 * Routine for parallel summation. Can be looped over from host.
 */
template<typename T> __global__ void pSumT(T* in1, T* output, int pass){
        unsigned int tid = threadIdx.x;
        unsigned int bid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x;// printf("bid0=%d\n",bid);
        unsigned int gid = getGid3d3d();
        extern __shared__ T sdata[];
        for(int i = blockDim.x>>1; i > 0; i>>=1){
                if(tid < blockDim.x>>1){
                        sdata[tid] += sdata[tid + i];
                }
                __syncthreads();
        }
        if(tid==0){
                output[bid] = sdata[0];
        }
}

/**
 * Routine for parallel summation. Can be looped over from host.
 */
__global__ void pSum(double* in1, double* output, int pass){
        unsigned int tid = threadIdx.x;
        unsigned int bid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x;// printf("bid0=%d\n",bid);
        unsigned int gid = getGid3d3d();
        extern __shared__ double sdata2[];
        for(int i = blockDim.x>>1; i > 0; i>>=1){
                if(tid < blockDim.x>>1){
                        sdata2[tid] += sdata2[tid + i];
                }
                __syncthreads();
        }
        if(tid==0){
                output[bid] = sdata2[0];
        }
}

//##############################################################################
//##############################################################################
