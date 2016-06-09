/*
* init.cu - GPUE: Split Operator based GPU solver for Nonlinear 
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
#include "../include/evolution.h"
#include <string>
#include <iostream>

/* Error variable & return variables */
cudaError_t err;
cufftResult result;

/* Define operating modes */
int ang_mom = 0;
int gpe = 0;

/* Allocating global variables */
double mass, a_s, omegaX, omegaY, omegaZ;
double xi; //Healing length minimum value defined at central density.

/* Evolution timestep */
double dt, gdt;

/* Grid dimensions vector. xyz are dim length, w is total grid size (x*y*z) */
int xDim, yDim, read_wfc, print, write_it;
long gsteps, esteps, atoms;
double *x,*y,*xp,*yp,*px,*py,dx,dy,xMax,yMax;

/* CuFFT plans for forward and inverse. May only need to use 1 for both */
cufftHandle plan_2d, plan_1d;

/* Arrays for storing wavefunction, momentum and position op, etc */
cufftDoubleComplex *wfc, *wfc0, *wfc_backup, *GK, *GV_half, *GV, *EK, *EV, *EV_opt, *GxPy, *GyPx, *ExPy, *EyPx, *EappliedField;
double *Energy, *Energy_gpu, *r, *Phi, *V, *V_opt, *K, *xPy, *yPx, *xPy_gpu, *yPx_gpu;

/* CUDA data buffers for FFT */
cufftDoubleComplex *wfc_gpu, *K_gpu, *V_gpu, *par_sum;
double *Phi_gpu;

/* CUDA streams */
cudaStream_t streamA, streamB, streamC, streamD;

/* Scaling the interaction */
double interaction;
double laser_power;

/* Define global dim3 and threads for grid and thread dim calculation */
dim3 grid;
int threads;

/* */
double l;

int initialise(Grid &par){

    char buffer[100];
    double gammaY; //Aspect ratio of trapping geometry.
    double Rxy; //Condensate scaling factor.
    double a0x, a0y; //Harmonic oscillator length in x and y directions

    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    unsigned int xD=1,yD=1,zD=1;
    threads = 128;

    // Re-establishing variables from parsed Grid class
    double omega = par.dval("omega");
    double angle_sweep = par.dval("angle_sweep");
    int kick_it = par.ival("kick_it");
    int graph = par.ival("graph");
    int N = par.ival("atoms");
    int printSteps = par.ival("print");
    int nonlin = par.ival("gpe");
    int lz = par.ival("ang_mom");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");

    // number of blocks in simulation
    unsigned int b = xDim*yDim/threads;

    // largest number of elements
    unsigned long long maxElements = 65536*65536ULL; 

    if( b < (1<<16) ){
        xD = b;
    }
    else if( (b >= (1<<16) ) && (b <= (maxElements)) ){
        int t1 = log(b)/log(2);
        float t2 = (float) t1/2;
        t1 = (int) t2;
        if(t2 > (float) t1){
            xD <<= t1;
            yD <<= (t1 + 1);
        }
        else if(t2 == (float) t1){
            xD <<= t1;
            yD <<= t1;
        }
    }
    else{
        printf("Outside range of supported indexing");
        exit(-1);
    }
    printf("Compute grid dimensions chosen as X=%d    Y=%d\n",xD,yD);
    
    grid.x=xD; 
    grid.y=yD; 
    grid.z=zD; 
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    int i,j; //Used in for-loops for indexing
    
    unsigned int gSize = xDim*yDim;
    double xOffset, yOffset;
    xOffset=0.0;//5.0e-6;
    yOffset=0.0;//5.0e-6;
    
    mass = 1.4431607e-25; //Rb 87 mass, kg
    par.store("Mass",mass);
    a_s = 4.67e-9;
    par.store("a_s",a_s);

    double sum = 0.0;

    a0x = sqrt(HBAR/(2*mass*omegaX));
    a0y = sqrt(HBAR/(2*mass*omegaY));
    par.store("a0x",a0x);
    par.store("a0y",a0y);
    
    Rxy = pow(15,0.2)*pow(N*a_s*sqrt(mass*omegaZ/HBAR),0.2);
    par.store("Rxy",Rxy);
    double bec_length = sqrt( HBAR/(mass*sqrt( omegaX*omegaX * 
                                               ( 1 - omega*omega) ) ));
    xMax = 6*Rxy*a0x; //10*bec_length; //6*Rxy*a0x;
    yMax = 6*Rxy*a0y; //10*bec_length;
    par.store("xMax",xMax);
    par.store("yMax",yMax);

    double pxMax, pyMax;
    pxMax = (PI/xMax)*(xDim>>1);
    pyMax = (PI/yMax)*(yDim>>1);
    par.store("pyMax",pyMax);
    par.store("pxMax",pxMax);
    
    dx = xMax/(xDim>>1);
    dy = yMax/(yDim>>1);
    par.store("dx",dx);
    par.store("dy",dy);
    
    double dpx, dpy;
    dpx = PI/(xMax);
    dpy = PI/(yMax);
    par.store("dpx",dpx);
    par.store("dpy",dpy);

    //printf("a0x=%e  a0y=%e \n dx=%e   dx=%e\n R_xy=%e\n",a0x,a0y,dx,dy,Rxy);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    //double *x,*y,*xp,*yp;
    x = (double *) malloc(sizeof(double) * xDim);
    y = (double *) malloc(sizeof(double) * yDim);
    xp = (double *) malloc(sizeof(double) * xDim);
    yp = (double *) malloc(sizeof(double) * yDim);

    /*
     * R-space and K-space grids
     */
    for(i=0; i<xDim/2; ++i){
        x[i] = -xMax + (i+1)*dx;        
        x[i + (xDim/2)] = (i+1)*dx;
        
        y[i] = -yMax + (i+1)*dy;        
        y[i + (yDim/2)] = (i+1)*dy;
        
        xp[i] = (i+1)*dpx;
        xp[i + (xDim/2)] = -pxMax + (i+1)*dpx;
        
        yp[i] = (i+1)*dpy;
        yp[i + (yDim/2)] = -pyMax + (i+1)*dpy;
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    /* Initialise wavefunction, momentum, position, angular momentum, 
       imaginary and real-time evolution operators . */
    Energy = (double*) malloc(sizeof(double) * gSize);
    r = (double *) malloc(sizeof(double) * gSize);
    Phi = (double *) malloc(sizeof(double) * gSize);
    wfc = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    wfc_backup = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * 
                                               (gSize/threads));
    K = (double *) malloc(sizeof(double) * gSize);
    V = (double *) malloc(sizeof(double) * gSize);
    V_opt = (double *) malloc(sizeof(double) * gSize);
    GK = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    GV = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EK = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EV = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EV_opt = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    xPy = (double *) malloc(sizeof(double) * gSize);
    yPx = (double *) malloc(sizeof(double) * gSize);
    ExPy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EyPx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EappliedField = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * 
                                                         gSize);
    
    /* Initialise wfc, EKp, and EVr buffers on GPU */
    cudaMalloc((void**) &Energy_gpu, sizeof(double) * gSize);
    cudaMalloc((void**) &wfc_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &Phi_gpu, sizeof(double) * gSize);
    cudaMalloc((void**) &K_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &V_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &xPy_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &yPx_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &par_sum, sizeof(cufftDoubleComplex) * (gSize/threads));
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    #ifdef __linux
    int cores = omp_get_num_procs();
    par.store("Cores_Total",cores);

    // Assuming dev system specifics (Xeon with HT -> cores detected / 2)
    par.store("Cores_Max",cores/2);
    omp_set_num_threads(cores/2);
    #pragma omp parallel for private(j)
    #endif
    for( i=0; i < xDim; i++ ){
        for( j=0; j < yDim; j++ ){
            Phi[(i*yDim + j)] = fmod(l*atan2(y[j], x[i]),2*PI);
            
            wfc[(i*yDim + j)].x = exp(-( pow((x[i])/(Rxy*a0x),2) + 
                                         pow((y[j])/(Rxy*a0y),2) ) ) *
                                  cos(Phi[(i*xDim + j)]);
            wfc[(i*yDim + j)].y = -exp(-( pow((x[i])/(Rxy*a0x),2) + 
                                          pow((y[j])/(Rxy*a0y),2) ) ) *
                                  sin(Phi[(i*xDim + j)]);
                
            V[(i*yDim + j)] = 0.5*mass*( pow(omegaX*(x[i]+xOffset),2) + 
                                         pow(gammaY*omegaY*(y[j]+yOffset),2) );
            K[(i*yDim + j)] = (HBAR*HBAR/(2*mass))*(xp[i]*xp[i] + yp[j]*yp[j]);

            GV[(i*yDim + j)].x = exp( -V[(i*xDim + j)]*(gdt/(2*HBAR)));
            GK[(i*yDim + j)].x = exp( -K[(i*xDim + j)]*(gdt/HBAR));
            GV[(i*yDim + j)].y = 0.0;
            GK[(i*yDim + j)].y = 0.0;
            
            xPy[(i*yDim + j)] = x[i]*yp[j];
            yPx[(i*yDim + j)] = -y[j]*xp[i];
            
            EV[(i*yDim + j)].x=cos( -V[(i*xDim + j)]*(dt/(2*HBAR)));
            EV[(i*yDim + j)].y=sin( -V[(i*xDim + j)]*(dt/(2*HBAR)));
            EK[(i*yDim + j)].x=cos( -K[(i*xDim + j)]*(dt/HBAR));
            EK[(i*yDim + j)].y=sin( -K[(i*xDim + j)]*(dt/HBAR));
            
            ExPy[(i*yDim + j)].x=cos(-omega*omegaX*xPy[(i*xDim + j)]*dt);
            ExPy[(i*yDim + j)].y=sin(-omega*omegaX*xPy[(i*xDim + j)]*dt);
            EyPx[(i*yDim + j)].x=cos(-omega*omegaX*yPx[(i*xDim + j)]*dt);
            EyPx[(i*yDim + j)].y=sin(-omega*omegaX*yPx[(i*xDim + j)]*dt);
    
            sum+=sqrt(wfc[(i*xDim + j)].x*wfc[(i*xDim + j)].x + 
                      wfc[(i*xDim + j)].y*wfc[(i*xDim + j)].y);
        }
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //hdfWriteDouble(xDim, V, 0, "V_0"); //HDF COMING SOON!
    //hdfWriteComplex(xDim, wfc, 0, "wfc_0");
    FileIO::writeOutDouble(buffer,"V",V,xDim*yDim,0);
    //FileIO::writeOutDouble(buffer,"V_opt",V_opt,xDim*yDim,0);
    FileIO::writeOutDouble(buffer,"K",K,xDim*yDim,0);
    FileIO::writeOutDouble(buffer,"xPy",xPy,xDim*yDim,0);
    FileIO::writeOutDouble(buffer,"yPx",yPx,xDim*yDim,0);
    FileIO::writeOut(buffer,"WFC",wfc,xDim*yDim,0);
    FileIO::writeOut(buffer,"ExPy",ExPy,xDim*yDim,0);
    FileIO::writeOut(buffer,"EyPx",EyPx,xDim*yDim,0);
    FileIO::writeOutDouble(buffer,"Phi",Phi,xDim*yDim,0);
    FileIO::writeOutDouble(buffer,"r",r,xDim*yDim,0);
    FileIO::writeOutDouble(buffer,"x",x,xDim,0);
    FileIO::writeOutDouble(buffer,"y",y,yDim,0);
    FileIO::writeOutDouble(buffer,"px",xp,xDim,0);
    FileIO::writeOutDouble(buffer,"py",yp,yDim,0);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //free(V); 
    free(K); free(r); //free(Phi);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    sum=sqrt(sum*dx*dy);
    //#pragma omp parallel for reduction(+:sum) private(j)
    for (i = 0; i < xDim; i++){
        for (j = 0; j < yDim; j++){
            wfc[(i*yDim + j)].x = (wfc[(i*yDim + j)].x)/(sum);
            wfc[(i*yDim + j)].y = (wfc[(i*yDim + j)].y)/(sum);
        }
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    result = cufftPlan2d(&plan_2d, xDim, yDim, CUFFT_Z2Z);
    if(result != CUFFT_SUCCESS){
        printf("Result:=%d\n",result);
        printf("Error: Could not execute cufftPlan2d(%s ,%d, %d).\n", "plan_2d",
                (unsigned int)xDim, (unsigned int)yDim);
        return -1;
    }

    result = cufftPlan1d(&plan_1d, xDim, CUFFT_Z2Z, yDim);
    if(result != CUFFT_SUCCESS){
        printf("Result:=%d\n",result);
        printf("Error: Could not execute cufftPlan3d(%s ,%d ,%d ).\n", 
               "plan_1d", (unsigned int)xDim, (unsigned int)yDim);
        return -1;
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    return 0;
}

// NOTE: RE-ESTABLISH PARAMS AFTER PARSING
int main(int argc, char **argv){

    char buffer[100];
    
    time_t start,fin;
    time(&start);
    printf("Start: %s\n", ctime(&start));
    //initArr(&params,32);
    //appendData(&params,ctime(&start),0.0);
    Grid par = parseArgs(argc,argv);
    Cuda cupar;
    int device = par.ival("device");
    cudaSetDevice(device);
    //************************************************************//
    /*
    * Initialise the Params data structure to track params and variables
    */
    //************************************************************//
    //paramS = (Params *) malloc(sizeof(Params));
    //strcpy(paramS->data,"INIT");
    //paramS->next=NULL;

    initialise(par);
    //************************************************************//
    /*
    * Groundstate finder section
    */
    //************************************************************//
    FileIO::writeOutParam(buffer, par, "Params.dat");
    if(read_wfc == 1){
        printf("Loading wavefunction...");
        wfc=FileIO::readIn("wfc_load","wfci_load",xDim, yDim);
        printf("Wavefunction loaded.\n");
    }
    
/*
    double x_0,y_0;
    x_0 = 0;//(0.5*xDim)*dx;
    y_0 = 0;//(0.5*yDim)*dy;
    for(int i=0; i < xDim; i++ ){
        for(int j=0; j < yDim; j++ ){
            ph.x = cos( fmod( 0*atan2( y[j] - y_0, x[i] - x_0 ), 2*PI) );
            ph.y = -sin( fmod( 0*atan2( y[j] - y_0, x[i] - x_0 ), 2*PI) );
            wfc[(i*yDim + j)] = Minions::complexMult( wfc[(i*yDim + j)], ph );
        }
    }
    printf("l=%e\n",l);
*/
    if(gsteps > 0){
        err=cudaMemcpy(K_gpu, GK, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(V_gpu, GV, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(xPy_gpu, xPy, sizeof(double)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(yPx_gpu, yPx, sizeof(double)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        
        evolve(wfc_gpu, K_gpu, V_gpu, yPx_gpu, xPy_gpu, par_sum, 
               par.ival("gsteps"), cupar, 0, 0, par, buffer);
        cudaMemcpy(wfc, wfc_gpu, sizeof(cufftDoubleComplex)*xDim*yDim, 
                   cudaMemcpyDeviceToHost);
    }

    free(GV); free(GK); free(xPy); free(yPx);

    //************************************************************//
    /*
    * Evolution
    */
    //************************************************************//
    if(esteps > 0){
        err=cudaMemcpy(xPy_gpu, ExPy, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(yPx_gpu, EyPx, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(xPy_gpu, ExPy, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(yPx_gpu, EyPx, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(K_gpu, EK, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim, 
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
            
        // delta_define(x, y, (523.6667 - 512 + x0_shift)*dx, 
        //              (512.6667 - 512 + y0_shift)*dy, V_opt);
        FileIO::writeOutDouble(buffer,"V_opt",V_opt,xDim*yDim,0);
        evolve(wfc_gpu, K_gpu, V_gpu, yPx_gpu, xPy_gpu, par_sum, 
               par.ival("esteps"), cupar, 1, 0, par, buffer);
    
    }
    free(EV); free(EK); free(ExPy); free(EyPx);
    free(x);free(y);
    cudaFree(wfc_gpu); cudaFree(K_gpu); cudaFree(V_gpu); cudaFree(yPx_gpu); 
    cudaFree(xPy_gpu); cudaFree(par_sum);

    time(&fin);
    //appendData(&params,ctime(&fin),0.0);
    printf("Finish: %s\n", ctime(&fin));
    printf("Total time: %ld seconds\n ",(long)fin-start);
    //appendData(&params,"t_duration",fin-start);
    return 0;
}
