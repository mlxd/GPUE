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

int initialise(Op &opr, Cuda &cupar, Grid &par, Wave &wave){

    // Re-establishing variables from parsed Grid class
    // Initializes uninitialized variables to 0 values
    int N = par.ival("atoms");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int threads = par.ival("threads");
    unsigned int gSize = xDim*yDim;
    double omega = par.dval("omega");
    double gdt = par.dval("gdt");
    double dt = par.dval("dt");
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double omegaZ = par.dval("omegaZ");
    double mass = par.dval("mass");
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double a_s = par.dval("a_s");
    double xMax = par.dval("xMax");
    double yMax = par.dval("ymax");
    double l = par.dval("l");
    double *x;
    double *y;
    double *xp;
    double *yp;
    double *Energy;
    double *r;
    double *V;
    double *V_opt;
    double *Phi;
    double *Phi_gpu;
    double *K;
    double *xPy;
    double *yPx;
    double *Energy_gpu;
    cufftDoubleComplex *wfc;
    //cufftDoubleComplex *V_gpu;
    cufftDoubleComplex *EV_opt;
    cufftDoubleComplex *wfc_backup;
    cufftDoubleComplex *GK;
    cufftDoubleComplex *GV;
    cufftDoubleComplex *EV;
    cufftDoubleComplex *EK;
    cufftDoubleComplex *ExPy;
    cufftDoubleComplex *EyPx;
    cufftDoubleComplex *EappliedField; 
    //cufftDoubleComplex *wfc_gpu;
    //cufftDoubleComplex *K_gpu;
    //cufftDoubleComplex *xPy_gpu;
    //cufftDoubleComplex *yPx_gpu;
    //cufftDoubleComplex *par_sum;

    cufftResult result = cupar.cufftResultval("result");
    cufftHandle plan_1d = cupar.cufftHandleval("plan_1d");
    cufftHandle plan_2d = cupar.cufftHandleval("plan_2d");

    dim3 grid = cupar.dim3val("grid");

    char buffer[100];
    double gammaY; //Aspect ratio of trapping geometry.
    double Rxy; //Condensate scaling factor.
    double a0x, a0y; //Harmonic oscillator length in x and y directions

    unsigned int xD=1,yD=1,zD=1;
    threads = 128;

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
    
    double xOffset, yOffset;
    xOffset=0.0;//5.0e-6;
    yOffset=0.0;//5.0e-6;
    
    mass = 1.4431607e-25; //Rb 87 mass, kg
    par.store("mass",mass);
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
    /*
    cudaMalloc((void**) &Energy_gpu, sizeof(double) * gSize);
    cudaMalloc((void**) &wfc_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &Phi_gpu, sizeof(double) * gSize);
    cudaMalloc((void**) &K_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &V_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &xPy_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &yPx_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &par_sum, sizeof(cufftDoubleComplex) * (gSize/threads));
    */
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
    //free(K); free(r); //free(Phi);

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

    // Storing variables that have been initialized
    // Re-establishing variables from parsed Grid class
    // Initializes uninitialized variables to 0 values
    par.store("omega", omega);
    par.store("gdt", gdt);
    par.store("dt", dt);
    par.store("omegaX", omegaX);
    par.store("omegaY", omegaY);
    par.store("omegaZ", omegaZ);
    par.store("dx", dx);
    par.store("dy", dy);
    par.store("xMax", xMax);
    par.store("ymax", yMax);
    par.store("l", l);
    par.store("x", x);
    par.store("y", y);
    par.store("xp", xp);
    par.store("yp", yp);
    wave.store("Energy", Energy);
    wave.store("r", r);
    opr.store("V", V);
    opr.store("V_opt", V_opt);
    wave.store("Phi", Phi);
    //wave.store("Phi_gpu", Phi_gpu);
    opr.store("K", K);
    opr.store("xPy", xPy);
    opr.store("yPx", yPx);
    //opr.store("Energy_gpu", Energy_gpu);
    par.store("atoms", N);
    par.store("xDim", xDim);
    par.store("yDim", yDim);
    par.store("threads", threads);
    wave.store("wfc", wfc);
    //opr.store("V_gpu", V_gpu);
    opr.store("EV_opt", EV_opt);
    wave.store("wfc_backup", wfc_backup);
    opr.store("GK", GK);
    opr.store("GV", GK);
    opr.store("EV", EV);
    opr.store("EK", EK);
    opr.store("ExPy", ExPy);
    opr.store("EyPx", EyPx);
    opr.store("EappliedField", EappliedField);
    //wave.store("wfc_gpu", wfc_gpu);
    //opr.store("K_gpu", K_gpu);
    //opr.store("xPy_gpu", xPy_gpu);
    //opr.store("yPx_gpu", yPx_gpu);
    //wave.store("par_sum", par_sum);

    cupar.store("result", result);
    cupar.store("plan_1d", plan_1d);
    cupar.store("plan_2d", plan_2d);

    cupar.store("grid", grid);

    return 0;
}

// NOTE: RE-ESTABLISH PARAMS AFTER PARSING
int main(int argc, char **argv){

    char buffer[100];

    Wave wave;
    Op opr;
    
    time_t start,fin;
    time(&start);
    printf("Start: %s\n", ctime(&start));
    Grid par = parseArgs(argc,argv);
    Cuda cupar;
    int device = par.ival("device");
    cudaSetDevice(device);

    //************************************************************//
    /*
    * Initialise the Params data structure to track params and variables
    */
    //************************************************************//

    initialise(opr, cupar, par, wave);

    // Re-establishing variables from parsed Grid class
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double *V_opt = opr.dsval("V_opt");
    double *xPy = opr.dsval("xPy");
    double *yPx = opr.dsval("yPx");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int read_wfc = par.ival("read_wfc");
    int gsteps = par.ival("gsteps");
    int esteps = par.ival("esteps");
    cufftDoubleComplex *wfc = wave.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *V_gpu = opr.cufftDoubleComplexval("V_gpu");
    cufftDoubleComplex *GK = opr.cufftDoubleComplexval("GK");
    cufftDoubleComplex *GV = opr.cufftDoubleComplexval("GV");
    cufftDoubleComplex *EV = opr.cufftDoubleComplexval("EV");
    cufftDoubleComplex *EK = opr.cufftDoubleComplexval("EK");
    cufftDoubleComplex *ExPy = opr.cufftDoubleComplexval("ExPy");
    cufftDoubleComplex *EyPx = opr.cufftDoubleComplexval("EyPx");
    cufftDoubleComplex *wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
    cufftDoubleComplex *K_gpu = opr.cufftDoubleComplexval("K_gpu");
    cufftDoubleComplex *xPy_gpu = opr.cufftDoubleComplexval("xPy_gpu");
    cufftDoubleComplex *yPx_gpu = opr.cufftDoubleComplexval("yPx_gpu");
    cufftDoubleComplex *par_sum = wave.cufftDoubleComplexval("par_sum");
    cudaError_t err = cupar.cudaError_tval("err");

    std::cout << "variables re-established" << '\n';
    std::cout << read_wfc << '\n';

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

    std::cout << "gsteps: " << gsteps << '\n';
    
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
        
        //evolve(wfc_gpu, K_gpu, V_gpu, yPx_gpu, xPy_gpu, par_sum, 
        //       par.ival("gsteps"), cupar, 0, 0, par, buffer);
        evolve(wave, opr, par_sum,
               par.ival("gsteps"), cupar, 0, 0, par, buffer);
        wfc = wave.cufftDoubleComplexval("wfc");
        wfc_gpu = wave.cufftDoubleComplexval("wfc");
        cudaMemcpy(wfc, wfc_gpu, sizeof(cufftDoubleComplex)*xDim*yDim, 
                   cudaMemcpyDeviceToHost);
    }

    std::cout << "got to here" << '\n';

    free(GV); free(GK); free(xPy); free(yPx);

    std::cout << "evolution started..." << '\n';
    std::cout << "esteps: " << esteps << '\n';

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
            
        FileIO::writeOutDouble(buffer,"V_opt",V_opt,xDim*yDim,0);
        //evolve(wfc_gpu, K_gpu, V_gpu, yPx_gpu, xPy_gpu, par_sum, 
        //       par.ival("esteps"), cupar, 1, 0, par, buffer);
        evolve(wave, opr, par_sum,
               par.ival("esteps"), cupar, 1, 0, par, buffer);
    
    }
    free(EV); free(EK); free(ExPy); free(EyPx);
    free(x);free(y);
    cudaFree(wfc_gpu); cudaFree(K_gpu); cudaFree(V_gpu); cudaFree(yPx_gpu); 
    cudaFree(xPy_gpu); cudaFree(par_sum);

    time(&fin);
    printf("Finish: %s\n", ctime(&fin));
    printf("Total time: %ld seconds\n ",(long)fin-start);
    return 0;
}
