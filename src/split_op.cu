/*
* split_op.cu - GPUE: Split Operator based GPU solver for Nonlinear 
Schrodinger Equation, Copyright (C) 2011-2015, Lee J. O'Riordan 
<loriordan@gmail.com>, Tadhg Morgan, Neil Crowley. All rights reserved.

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

#include "../include/split_op.h"
#include "../include/kernels.h"
#include "../include/constants.h"
#include "../include/fileIO.h"
#include "../include/tracker.h"
#include "../include/minions.h"
#include "../include/parser.h"

#include "../include/lattice.h"
#include "../include/node.h"
#include "../include/edge.h"
#include "../include/manip.h"
#include "../include/vort.h"
#include <string>
#include <iostream>

unsigned int LatticeGraph::Edge::suid = 0;
unsigned int LatticeGraph::Node::suid = 0;

char buffer[100];
int verbose; //Print more info. Not curently implemented.
int device; //GPU ID choice.
int kick_it; //Kicking mode: 0 = off, 1 = multiple, 2 = single
int graph=0; //Generate graph from vortex lattice.
double gammaY; //Aspect ratio of trapping geometry.
double omega; //Rotation rate of condensate
double timeTotal;
double angle_sweep; //Rotation angle of condensate relative to x-axis
double x0_shift, y0_shift; //Optical lattice shift parameters.
double Rxy; //Condensate scaling factor.
double a0x, a0y; //Harmonic oscillator length in x and y directions
double sepMinEpsilon=0.0; //Minimum separation for epsilon.

/*
 * Checks CUDA routines have exitted correctly.
 */
int isError(int result, char* c){
    if(result!=0){
        printf("Error has occurred for method %s with return type %d\n",
               c,result);
        exit(result);
    }
    return result;
}

/*
 * Used to perform parallel summation on WFC for normalisation.
 */
void parSum(double2* gpuWfc, double2* gpuParSum, Grid &par, 
            Cuda &cupar){ 
    // May need to add double l
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    int threads = par.ival("threads");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int grid_tmp = xDim*yDim;
    int block = grid_tmp/threads;
    int thread_tmp = threads;
    int pass = 0;
    dim3 grid = cupar.dim3val("grid");
    while((double)grid_tmp/threads > 1.0){
        if(grid_tmp == xDim*yDim){
            multipass<<<block,threads,threads*sizeof(double2)>>>(&gpuWfc[0],
                &gpuParSum[0],pass); 
        }
        else{
            multipass<<<block,thread_tmp,thread_tmp*sizeof(double2)>>>(
                &gpuParSum[0],&gpuParSum[0],pass);
        }
        grid_tmp /= threads;
        block = (int) ceil((double)grid_tmp/threads);
        pass++;
    }
    thread_tmp = grid_tmp;
    multipass<<<1,thread_tmp,thread_tmp*sizeof(double2)>>>(&gpuParSum[0],
                                                           &gpuParSum[0], pass);
    scalarDiv_wfcNorm<<<grid,threads>>>(gpuWfc, dx*dy, gpuParSum, gpuWfc);
}

/**
** Matches the optical lattice to the vortex lattice. 
** Moire super-lattice project.
**/
void optLatSetup(struct Vtx::Vortex centre, double* V, 
                 struct Vtx::Vortex *vArray, int num_vortices, double theta_opt,
                 double intensity, double* v_opt, double *x, double *y,
                 Grid &par, Op &opr){
    std::string data_dir = par.sval("data_dir");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double dt = par.dval("dt");
    cufftDoubleComplex *EV_opt = opr.cufftDoubleComplexval("EV_opt");
    int i,j;
    double sepMin = Tracker::vortSepAvg(vArray,centre,num_vortices);
    sepMin = sepMin*(1 + sepMinEpsilon);
    par.store("Vort_sep",(double)sepMin);
    
    // Defining the necessary k vectors for the optical lattice
    

    // Additional /2 as a result of lambda/2 period
    double k_mag = ((2*PI/(sepMin*dx))/2)*(2/sqrt(3));
    double2* k = (double2*) malloc(sizeof(double2)*3);
    par.store("kmag",(double)k_mag);
    k[0].x = k_mag * cos(0*PI/3 + theta_opt);
    k[0].y = k_mag * sin(0*PI/3 + theta_opt);
    k[1].x = k_mag * cos(2*PI/3 + theta_opt);
    k[1].y = k_mag * sin(2*PI/3 + theta_opt);
    k[2].x = k_mag * cos(4*PI/3 + theta_opt);
    k[2].y = k_mag * sin(4*PI/3 + theta_opt);
    
    double2 *r_opt = (double2*) malloc(sizeof(double2)*xDim);

    //FileIO::writeOut(buffer,data_dir + "r_opt",r_opt,xDim,0);
    par.store("k[0].x",(double)k[0].x);
    par.store("k[0].y",(double)k[0].y);
    par.store("k[1].x",(double)k[1].x);
    par.store("k[1].y",(double)k[1].y);
    par.store("k[2].x",(double)k[2].x);
    par.store("k[2].y",(double)k[2].y);

    // sin(theta_opt)*(sepMin);
    double x_shift = dx*(9+(0.5*xDim-1) - centre.coords.x);

    // cos(theta_opt)*(sepMin);
    double y_shift = dy*(0+(0.5*yDim-1) - centre.coords.y);

    printf("Xs=%e\nYs=%e\n",x_shift,y_shift);

    //#pragma omp parallel for private(j)
    for ( j=0; j<yDim; ++j ){
        for ( i=0; i<xDim; ++i ){
            v_opt[j*xDim + i] = intensity*(
                                pow( ( cos( k[0].x*( x[i] + x_shift ) + 
                                       k[0].y*( y[j] + y_shift ) ) ), 2) +
                                pow( ( cos( k[1].x*( x[i] + x_shift ) + 
                                       k[1].y*( y[j] + y_shift ) ) ), 2) +
                                pow( ( cos( k[2].x*( x[i] + x_shift ) + 
                                       k[2].y*( y[j] + y_shift ) ) ), 2)
                                );
            EV_opt[(j*xDim + i)].x=cos( -(V[(j*xDim + i)] + 
                                   v_opt[j*xDim + i])*(dt/(2*HBAR)));
            EV_opt[(j*xDim + i)].y=sin( -(V[(j*xDim + i)] + 
                                   v_opt[j*xDim + i])*(dt/(2*HBAR)));
        }
    }

    // Storing changed variables
    opr.store("EV_opt", EV_opt);
    opr.store("V", V);
    opr.store("V_opt",v_opt);
}

/**
** Calculates energy and angular momentum of current state. 
** Implementation not fully finished.
**/
/*
double energy_angmom(double *Energy, double* Energy_gpu, double2 *V_op, 
                     double2 *K_op, double2 *gpuWfc, 
                     int gState, Grid &par){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double dt = par.dval("dt");

    double renorm_factor_2d=1.0/pow(xDim*yDim,0.5);
    double result=0;

    for (int i=0; i < xDim*yDim; ++i){
        Energy[i] = 0.0; 
    }
    return result*dx*dy;
    
}

// Creates narrow Gaussian "delta" peaks for vortex kicking
void delta_define(double *x, double *y, double x0, double y0, double *delta,
                  Grid &par, Op &opr){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    cufftDoubleComplex *EV_opt = opr.cufftDoubleComplexval("EV_opt");
    double *V = opr.dsval("V");
    double dx = par.dval("dx");
    double dt = par.dval("dt");

    for (int i=0; i<xDim; ++i){
        for (int j=0; j<yDim; ++j){
            delta[j*xDim + i] = 1e6*HBAR*exp( -( pow( x[i] - x0, 2) + 
                                pow( y[j] - y0, 2) )/(5*dx*dx) );
            EV_opt[(j*xDim + i)].x = cos( -(V[(j*xDim + i)] + 
                                     delta[j*xDim + i])*(dt/(2*HBAR)));
            EV_opt[(j*xDim + i)].y = sin( -(V[(j*xDim + i)] + 
                                     delta[j*xDim + i])*(dt/(2*HBAR)));
        }
    }
}
*/
