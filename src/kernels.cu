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


__constant__ double gDenConst = 6.6741e-40;//Evaluted in MATLAB: N*4*HBAR*HBAR*PI*(4.67e-9/mass)*sqrt(mass*(omegaZ)/(2*PI*HBAR))
//inline __device__ unsigned int getGid3d3d(){

__device__ unsigned int getGid3d3d(){
	return blockDim.x * ( ( blockDim.y * ( ( blockIdx.z * blockDim.z + threadIdx.z ) + blockIdx.y ) + threadIdx.y ) + blockIdx.x ) + threadIdx.x;
}

//inline __device__ unsigned int getBid3d3d(){
__device__ unsigned int getBid3d3d(){
	return blockIdx.x + gridDim.x*(blockIdx.y + gridDim.y * blockIdx.z);
}


//inline __device__ unsigned int getTid3d3d(){
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

//inline __device__ double complexMagnitude(double2 in){
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
	double2 result;
	unsigned int gid = getGid3d3d();
	result.x = (in1[gid].x*in2[gid].x - in1[gid].y*in2[gid].y);
	result.y = (in1[gid].x*in2[gid].y + in1[gid].y*in2[gid].x);
	out[gid] = result;
}

__global__ void cMultDensity(double2* in1, double2* in2, double2* out, double dt, double mass,double omegaZ, int gstate, int N){
	double2 result;
	double gDensity;
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	gDensity = gDenConst*complexMagnitudeSquared(in2[tid]); // scaling of interaction strength doesn't work now

	if(gstate == 0){
		double tmp = in1[tid].x*exp(-gDensity*(dt/HBAR) );
		result.x = (tmp)*in2[tid].x - (in1[tid].y)*in2[tid].y;
		result.y = (tmp)*in2[tid].y + (in1[tid].y)*in2[tid].x;
	}
	else{
		double2 tmp;
		tmp.x = in1[tid].x*cos(-gDensity*(dt/HBAR)) - in1[tid].y*sin(-gDensity*(dt/HBAR));
		tmp.y = in1[tid].y*cos(-gDensity*(dt/HBAR)) + in1[tid].x*sin(-gDensity*(dt/HBAR));
		
		result.x = (tmp.x)*in2[tid].x - (tmp.y)*in2[tid].y;
		result.y = (tmp.x)*in2[tid].y + (tmp.y)*in2[tid].x;
	}
	out[tid] = result;
}

/**
 * Divides both components of vector type "in", by the value "factor".
 * Results given with "out"
 */
__global__ void scalarDiv(double2* in, double factor, double2* out){
	double2 result;
	//extern __shared__ double2 tmp_in[];
	unsigned int gid = getGid3d3d();
	result.x = (in[gid].x*factor);
	result.y = (in[gid].y*factor);
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
 */
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
	unsigned int bid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x;// printf("bid0=%d\n",bid);
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
	//printf("oplocal=%e	Resx=%e	Resy=%e\n",opLocal,result.x,result.y);
	energy[gid].x += result.x;
	energy[gid].y += result.y;
}


//#############################################################################################
//#############################################################################################

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



//#############################################################################################
//#############################################################################################
