/*
* kernels.cu - GPUE: Split Operator based GPU solver for Nonlinear 
* Schrodinger Equation, Copyright (C) 2012, Lee J. O'Riordan, Tadhg 
* Morgan, Neil Crowley. 

* This library is free software; you can redistribute it and/or modify 
* it under the terms of the GNU Lesser General Public License as 
* published by the Free Software Foundation; either version 2.1 of the 
* License, or (at your option) any later version. This library is 
* distributed in the hope that it will be useful, but WITHOUT ANY 
* WARRANTY; without even the implied warranty of MERCHANTABILITY or 
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
* License for more details. You should have received a copy of the GNU 
* Lesser General Public License along with this library; if not, write 
* to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
* Boston, MA 02111-1307 USA 
*/

#include "../include/constants.h"

/* Constant memory */
__const__ double FFT_RENORM_2D = 0.0;
__const__ double FFT_RENORM_1D = 0.0;

__device__ double complexMagnitude(double2 in){
	return sqrt(in.x*in.x + in.y*in.y);
}

__device__ double complexMagnitudeSquared(double2 in){
	return in.x*in.x + in.y*in.y;
}

__device__ double2 complexMultiply(double2 in1, double2 in2){
	double2 result;
	result.x = (in1.x*in2.x - in1.y*in2.y);
    result.y = (in1.x*in2.y + in1.y*in2.x);
	return result;
}

/**
 * Used to pin a vortex to a specified location on the condensate
 */
__global__ void pinVortex(double2* wfc_in, double2* wfc_out, double2* r){
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	double radius = sqrt( (r[tid].x)*(r[tid].x) + (r[tid].y)*(r[tid].y) );
	double2 multFactor;
	if(radius == 0.0){
		multFactor.x = 1.0;
		multFactor.y = 0.0;
	}
	else{
		multFactor.x = (r[tid].x)/radius;
		multFactor.y = (r[tid].y)/radius;
	}
	wfc_out[tid] = complexMultiply(wfc_in[tid], multFactor);
}


/**
 * Naive transpose in global memory only
 */
__global__ void transpose(double2 *in, double2* out){
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	int yIdx = blockIdx.y*gridDim.x*blockDim.x;
	int xIdx = blockIdx.x*blockDim.x + threadIdx.x;
	//__shared__ double2 mat1[][];
	//__shared__ double2 mat2[][];
}

/**
 * Performs complex multiplication of in1 and in2, giving result as out. 
 */
/*
__global__ void cMult(double2* in1, double2* in2, double2* out){
	double2 result;
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	extern __shared__ double2 Mat1[];
	extern __shared__ double2 Mat2[];
	Mat1[threadIdx.x] = in1[tid];
	Mat2[threadIdx.x] = in2[tid];
	__syncthreads();
	result.x = (Mat1[threadIdx.x].x*Mat2[threadIdx.x].x - Mat1[threadIdx.x].y*Mat2[threadIdx.x].y);
	result.y = (Mat1[threadIdx.x].x*Mat2[threadIdx.x].y + Mat1[threadIdx.x].y*Mat2[threadIdx.x].x);
	__syncthreads();
	out[tid] = result;
}*/


/**
 * Performs complex multiplication of in1 and in2, giving result as out. 
 */
__global__ void cMult(double2* in1, double2* in2, double2* out){
	double2 result;
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	result.x = (in1[tid].x*in2[tid].x - in1[tid].y*in2[tid].y);
	result.y = (in1[tid].x*in2[tid].y + in1[tid].y*in2[tid].x);
	out[tid] = result;
}

__global__ void cMultDensity(double2* in1, double2* in2, double2* out, double dt, double mass,double omegaZ, int gstate, int N){
	double2 result;
	double gDensity;
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	gDensity = (0.5*N)*complexMagnitudeSquared(in2[tid])*4*HBAR*HBAR*PI*(4.67e-9/mass)*sqrt(mass*(omegaZ)/(2*PI*HBAR));

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
__global__ void scalarDiv2D(double2* in, double2* out){
	double2 result;
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	result.x = (in[tid].x*FFT_RENORM_2D);
	result.y = (in[tid].y*FFT_RENORM_2D);
	out[tid] = result;
}

/**
 * Divides both components of vector type "in", by the value "factor".
 * Results given with "out"
 */
__global__ void scalarDiv(double2* in, double factor, double2* out){
	double2 result;
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	result.x = (in[tid].x*factor);
	result.y = (in[tid].y*factor);
	out[tid] = result;
}

__global__ void scalarDiv_wfcNorm(double2* in, double dr, double2* pSum, double2* out){
	double2 result;
	int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	double norm = sqrt((pSum[0].x + pSum[0].y)*dr);
	result.x = (in[tid].x/norm);
	result.y = (in[tid].y/norm);
	out[tid] = result;
}

/**
 * Routine for parallel summation. Can be looped over from host.
 */
__global__ void multipass(double2* input, double2* output, int pass){
	unsigned int tid = threadIdx.x;
	unsigned int bid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x;
	unsigned int gid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
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
	}/*
	if(tid < 32){
		sdata[tid].x += sdata[tid + 32].x;
		sdata[tid].x += sdata[tid + 16].x;
		sdata[tid].x += sdata[tid + 8].x;
		sdata[tid].x += sdata[tid + 4].x;
		sdata[tid].x += sdata[tid + 2].x;
		sdata[tid].x += sdata[tid + 1].x;
		
		sdata[tid].y += sdata[tid + 32].y;
		sdata[tid].y += sdata[tid + 16].y;
		sdata[tid].y += sdata[tid + 8].y;
		sdata[tid].y += sdata[tid + 4].y;
		sdata[tid].y += sdata[tid + 2].y;
		sdata[tid].y += sdata[tid + 1].y;
	}*/
	if(tid==0){
		output[bid] = sdata[0];
	}
}

/*
__global__ void gpeEnergy(double2* V, double2* K, double2* wfc, double dr){
	uint tid = threadIdx.x;
	uint bid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x;
	uint gid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	
	//double E_K = 
}
*/
