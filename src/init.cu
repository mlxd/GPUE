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

#include "../include/init.h"

int init_2d(Op &opr, Cuda &cupar, Grid &par, Wave &wave){

    // Setting functions for operators
    opr.set_fns();

    // Re-establishing variables from parsed Grid class
    // Initializes uninitialized variables to 0 values
    std::string data_dir = par.sval("data_dir");
    int N = par.ival("atoms");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int threads;
    unsigned int gSize = xDim*yDim;
    double omega = par.dval("omega");
    double gdt = par.dval("gdt");
    double dt = par.dval("dt");
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double omegaZ = par.dval("omegaZ");
    double gammaY = par.dval("gammaY"); //Aspect ratio of trapping geometry.
    double l = par.dval("winding");
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
    double *pAy;
    double *pAx;
    double *Ax;
    double *Ay;
    double *Bz;
    double *pAy_gpu;
    double *pAx_gpu;
    double *Energy_gpu;
    cufftDoubleComplex *wfc;
    cufftDoubleComplex *V_gpu;
    cufftDoubleComplex *EV_opt;
    cufftDoubleComplex *wfc_backup;
    cufftDoubleComplex *GK;
    cufftDoubleComplex *GV;
    cufftDoubleComplex *GpAx;
    cufftDoubleComplex *GpAy;
    cufftDoubleComplex *EV;
    cufftDoubleComplex *EK;
    cufftDoubleComplex *EpAy;
    cufftDoubleComplex *EpAx;
    cufftDoubleComplex *EappliedField; 
    cufftDoubleComplex *wfc_gpu;
    cufftDoubleComplex *K_gpu;
    cufftDoubleComplex *par_sum;

    //std::cout << omegaX << '\t' << omegaY << '\n';
    //std::cout << "xDim is: " << xDim << '\t' <<  "yDim is: " << yDim << '\n';

    cufftResult result = cupar.cufftResultval("result");
    cufftHandle plan_1d;
    cufftHandle plan_2d;
    cufftHandle plan_other2d;

    dim3 grid = cupar.dim3val("grid");

    std::string buffer;
    double Rxy; //Condensate scaling factor.
    double a0x, a0y; //Harmonic oscillator length in x and y directions

    unsigned int xD=1,yD=1;
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
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    int i,j; //Used in for-loops for indexing
    
/*
    double xOffset, yOffset;
    xOffset=0.0;//5.0e-6;
    yOffset=0.0;//5.0e-6;
*/
    
    double mass = 1.4431607e-25; //Rb 87 mass, kg
    par.store("mass",mass);
    double a_s = 4.76e-9;
    par.store("a_s",a_s);

    double sum = 0.0;

    a0x = sqrt(HBAR/(2*mass*omegaX));
    a0y = sqrt(HBAR/(2*mass*omegaY));
    par.store("a0x",a0x);
    par.store("a0y",a0y);

    //std::cout << "a0x and y are: " << a0x << '\t' << a0y << '\n';

    //std::cout << N << '\t' << a_s << '\t' << mass << '\t' << omegaZ << '\n';
    
    Rxy = pow(15,0.2)*pow(N*a_s*sqrt(mass*omegaZ/HBAR),0.2);
    par.store("Rxy",Rxy);
    double bec_length = sqrt( HBAR/(mass*sqrt( omegaX*omegaX * 
                                               ( 1 - omega*omega) ) ));

    //std::cout << "Rxy is: " << Rxy << '\n';
    double xMax = 6*Rxy*a0x; //10*bec_length; //6*Rxy*a0x;
    double yMax = 6*Rxy*a0y; //10*bec_length;
    par.store("xMax",xMax);
    par.store("yMax",yMax);

    double pxMax, pyMax;
    pxMax = (PI/xMax)*(xDim>>1);
    pyMax = (PI/yMax)*(yDim>>1);
    par.store("pyMax",pyMax);
    par.store("pxMax",pxMax);
    
    double dx = xMax/(xDim>>1);
    double dy = yMax/(yDim>>1);
    par.store("dx",dx);
    par.store("dy",dy);
    
    double dpx, dpy;
    dpx = PI/(xMax);
    dpy = PI/(yMax);
    //std::cout << "yMax is: " << yMax << '\t' << "xMax is: " << xMax << '\n';
    //std::cout << "dpx and dpy are:" << '\n';
    //std::cout << dpx << '\t' << dpy << '\n';
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
    //std::cout << "dx and dy are: " << '\n';
    //std::cout << dx << '\t' << dy << '\n';

    // creating x,y,xp,yp
    for(i=0; i<xDim/2; ++i){
        x[i] = -xMax + i*dx;        
        x[i + (xDim/2)] = i*dx;
        
        xp[i] = i*dpx;
        xp[i + (xDim/2)] = -pxMax + i*dpx;
        
    }
    for(i=0; i<yDim/2; ++i){
        y[i] = -yMax + i*dy;        
        y[i + (yDim/2)] = i*dy;
        
        yp[i] = i*dpy;
        yp[i + (yDim/2)] = -pyMax + i*dpy;

    }

    par.store("x", x);
    par.store("y", y);
    par.store("xp", xp);
    par.store("yp", yp);
    

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
    GpAx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    GpAy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EK = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EV = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EV_opt = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    pAy = (double *) malloc(sizeof(double) * gSize);
    pAx = (double *) malloc(sizeof(double) * gSize);
    Ax = (double *) malloc(sizeof(double) * gSize);
    Ay = (double *) malloc(sizeof(double) * gSize);
    Bz = (double *) malloc(sizeof(double) * gSize);
    EpAy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EpAx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EappliedField = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * 
                                                         gSize);
    
    /* Initialise wfc, EKp, and EVr buffers on GPU */
    cudaMalloc((void**) &Energy_gpu, sizeof(double) * gSize);
    cudaMalloc((void**) &wfc_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &Phi_gpu, sizeof(double) * gSize);
    cudaMalloc((void**) &K_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &V_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &pAy_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &pAx_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &par_sum, sizeof(cufftDoubleComplex) * (gSize/threads));
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //std::cout << "all variables malloc'd" << '\n';

    #ifdef __linux
    int cores = omp_get_num_procs();
    par.store("Cores_Total",cores);

    // Assuming dev system specifics (Xeon with HT -> cores detected / 2)
    par.store("Cores_Max",cores/2);
    omp_set_num_threads(cores/2);

    // Setting Ax, and Ay if from file
    if (par.Afn == "file"){
        file_A(par.Axfile, Ax);
        opr.store("Ax",Ax);
        file_A(par.Ayfile, Ay);
        opr.store("Ay", Ay);
    }
    std::cout << "finished reading Ax / Ay from file" << '\n';
    std::cout << Ax[256] << '\t' << Ay[0] << '\n';
    #pragma omp parallel for private(j)
    #endif
    for( i=0; i < xDim; i++ ){
        for( j=0; j < yDim; j++ ){
            Phi[(i*yDim + j)] = fmod(l*atan2(y[j], x[i]),2*PI);
            
            if (par.bval("unit_test")){
                wfc[(i*yDim + j)].x =  (1/sqrt(2))*pow(1/PI,0.5) 
                    * exp( -0.5*( x[i]*x[i] + y[j]*y[j] ) )*(1+2*x[i]/sqrt(2));
                wfc[(i*yDim + j)].y = 0;
            }
            else if (par.bval("dimensionless")){
                wfc[(i*yDim + j)].x = exp(-( pow((x[i]),2) + 
                                             pow((y[j]),2) ) ) *
                                      cos(Phi[(i*xDim + j)]);
                wfc[(i*yDim + j)].y = -exp(-( pow((x[i]),2) + 
                                              pow((y[j]),2) ) ) *
                                          sin(Phi[(i*xDim + j)]);
            }
            else{
                wfc[(i*yDim + j)].x = exp(-( pow((x[i])/(Rxy*a0x),2) + 
                                             pow((y[j])/(Rxy*a0y),2) ) ) *
                                      cos(Phi[(i*xDim + j)]);
                wfc[(i*yDim + j)].y = -exp(-( pow((x[i])/(Rxy*a0x),2) + 
                                              pow((y[j])/(Rxy*a0y),2) ) ) *
                                          sin(Phi[(i*xDim + j)]);
            }
                
            V[(i*yDim + j)] = opr.V_fn(par.Vfn)(par, opr, i, j, 0);
            K[(i*yDim + j)] = opr.K_fn(par.Kfn)(par, opr, i, j, 0);

            GV[(i*yDim + j)].x = exp( -V[(i*xDim + j)]*(gdt/(2*HBAR)));
            GK[(i*yDim + j)].x = exp( -K[(i*xDim + j)]*(gdt/HBAR));
            GV[(i*yDim + j)].y = 0.0;
            GK[(i*yDim + j)].y = 0.0;

            // Ax and Ay will be calculated here but are used only for
            // debugging. They may be needed later for magnetic field calc
            if (par.Afn != "file"){
                Ax[(i*yDim + j)] = opr.Ax_fn(par.Afn)(par, opr, i, j, 0);
                Ay[(i*yDim + j)] = opr.Ay_fn(par.Afn)(par, opr, i, j, 0);
            }
            
            //pAy[(i*yDim + j)] = x[i]*yp[j];
            pAy[(i*yDim + j)] = opr.pAy_fn("rotation")(par, opr, i, j, 0);
            //pAx[(i*yDim + j)] = -y[j]*xp[i];
            pAx[(i*yDim + j)] = opr.pAx_fn("rotation")(par, opr, i, j, 0);

            GpAx[(i*yDim + j)].x = exp(-pAx[(i*xDim + j)]*gdt);
            GpAx[(i*yDim + j)].y = 0;
            GpAy[(i*yDim + j)].x = exp(-pAy[(i*xDim + j)]*gdt);
            GpAy[(i*yDim + j)].y = 0;
            
            EV[(i*yDim + j)].x=cos( -V[(i*xDim + j)]*(dt/(2*HBAR)));
            EV[(i*yDim + j)].y=sin( -V[(i*xDim + j)]*(dt/(2*HBAR)));
            EK[(i*yDim + j)].x=cos( -K[(i*xDim + j)]*(dt/HBAR));
            EK[(i*yDim + j)].y=sin( -K[(i*xDim + j)]*(dt/HBAR));
            
            EpAy[(i*yDim + j)].x=cos(-omega*omegaX*pAy[(i*xDim + j)]*dt);
            EpAy[(i*yDim + j)].y=sin(-omega*omegaX*pAy[(i*xDim + j)]*dt);
            EpAx[(i*yDim + j)].x=cos(-omega*omegaX*pAx[(i*xDim + j)]*dt);
            EpAx[(i*yDim + j)].y=sin(-omega*omegaX*pAx[(i*xDim + j)]*dt);
    
            sum+=sqrt(wfc[(i*xDim + j)].x*wfc[(i*xDim + j)].x + 
                      wfc[(i*xDim + j)].y*wfc[(i*xDim + j)].y);
        }
    }

    Bz = curl2d(par, Ax, Ay);

    std::cout << "writing initial variables to file..." << '\n';
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //hdfWriteDouble(xDim, V, 0, "V_0"); //HDF COMING SOON!
    //hdfWriteComplex(xDim, wfc, 0, "wfc_0");
    //FileIO::writeOutDouble(buffer, data_dir + "V_opt",V_opt,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "V",V,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "K",K,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "pAy",pAy,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "pAx",pAx,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "Ax",Ax,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "Ay",Ay,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "Bz",Bz,xDim*yDim,0);
    FileIO::writeOut(buffer, data_dir + "WFC",wfc,xDim*yDim,0);
    FileIO::writeOut(buffer, data_dir + "EpAy",EpAy,xDim*yDim,0);
    FileIO::writeOut(buffer, data_dir + "EpAx",EpAx,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "Phi",Phi,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "r",r,xDim*yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "x",x,xDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "y",y,yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "px",xp,xDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "py",yp,yDim,0);
    FileIO::writeOut(buffer, data_dir + "GK",GK,xDim*yDim,0);
    FileIO::writeOut(buffer, data_dir + "GV",GV,xDim*yDim,0);
    FileIO::writeOut(buffer, data_dir + "GpAx",GpAx,xDim*yDim,0);
    FileIO::writeOut(buffer, data_dir + "GpAy",GpAy,xDim*yDim,0);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //std::cout << "wrote initial variables" << '\n';

    //free(V); 
    free(K); free(r); free(Phi);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    sum=sqrt(sum*dx*dy);
    //#pragma omp parallel for reduction(+:sum) private(j)
    for (i = 0; i < xDim; i++){
        for (j = 0; j < yDim; j++){
            wfc[(i*yDim + j)].x = (wfc[(i*yDim + j)].x)/(sum);
            wfc[(i*yDim + j)].y = (wfc[(i*yDim + j)].y)/(sum);
        }
    }
    
    //std::cout << "modified wfc" << '\n';
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    //std::cout << "xDim is: " << xDim << '\t' << "yDim is: " << yDim << '\n';
    //std::cout << "plan_2d is: " << plan_2d << '\n';
    result = cufftPlan2d(&plan_2d, xDim, yDim, CUFFT_Z2Z);
    //std::cout << "found result" << '\n';
    if(result != CUFFT_SUCCESS){
        printf("Result:=%d\n",result);
        printf("Error: Could not execute cufftPlan2d(%s ,%d, %d).\n", "plan_2d",
                (unsigned int)xDim, (unsigned int)yDim);
        return -1;
    }

    plan_other2d = generate_plan_other2d(par); 

    result = cufftPlan1d(&plan_1d, xDim, CUFFT_Z2Z, yDim);
    if(result != CUFFT_SUCCESS){
        printf("Result:=%d\n",result);
        printf("Error: Could not execute cufftPlan1d(%s ,%d ,%d ).\n", 
               "plan_1d", (unsigned int)xDim, (unsigned int)yDim);
        return -1;
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //std::cout << GV[0].x << '\t' << GK[0].x << '\t' 
    //          << pAy[0] << '\t' << pAx[0] << '\n';

    //std::cout << "storing variables..." << '\n';

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
    par.store("yMax", yMax);
    par.store("winding", l);
    par.store("x", x);
    par.store("y", y);
    par.store("xp", xp);
    par.store("yp", yp);
    wave.store("Energy", Energy);
    wave.store("r", r);
    opr.store("V", V);
    opr.store("V_opt", V_opt);
    wave.store("Phi", Phi);
    wave.store("Phi_gpu", Phi_gpu);
    opr.store("K", K);
    opr.store("pAy", pAy);
    opr.store("pAx", pAx);
    opr.store("Energy_gpu", Energy_gpu);
    par.store("atoms", N);
    par.store("xDim", xDim);
    par.store("yDim", yDim);
    par.store("threads", threads);
    wave.store("wfc", wfc);
    opr.store("V_gpu", V_gpu);
    opr.store("EV_opt", EV_opt);
    wave.store("wfc_backup", wfc_backup);
    opr.store("GK", GK);
    opr.store("GV", GV);
    opr.store("GpAx", GpAx);
    opr.store("GpAy", GpAy);
    opr.store("EV", EV);
    opr.store("EK", EK);
    opr.store("EpAy", EpAy);
    opr.store("EpAx", EpAx);
    opr.store("EappliedField", EappliedField);
    wave.store("wfc_gpu", wfc_gpu);
    opr.store("K_gpu", K_gpu);
    opr.store("pAy_gpu", pAy_gpu);
    opr.store("pAx_gpu", pAx_gpu);
    wave.store("par_sum", par_sum);

    cupar.store("result", result);
    cupar.store("plan_1d", plan_1d);
    cupar.store("plan_2d", plan_2d);
    cupar.store("plan_other2d", plan_other2d);

    cupar.store("grid", grid);

    std::cout << "variables stored" << '\n';

    return 0;
}

// initializing all variables for 3d
int init_3d(Op &opr, Cuda &cupar, Grid &par, Wave &wave){

    // Setting functions for operators
    opr.set_fns();

    // Re-establishing variables from parsed Grid class
    // Initializes uninitialized variables to 0 values
    std::string data_dir = par.sval("data_dir");
    int N = par.ival("atoms");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int threads;
    unsigned int gSize = xDim*yDim*zDim;
    double omega = par.dval("omega");
    double gdt = par.dval("gdt");
    double dt = par.dval("dt");
    double omegaX = par.dval("omegaX");
    double omegaY = par.dval("omegaY");
    double omegaZ = par.dval("omegaZ");
    double gammaY = par.dval("gammaY"); //Aspect ratio of trapping geometry.
    double l = par.dval("winding");
    double *x;
    double *y;
    double *z;
    double *xp;
    double *yp;
    double *zp;
    double *Energy;
    double *r;
    double *V;
    double *V_opt;
    double *Phi;
    double *Phi_gpu;
    double *K;
    double *pAy;
    double *pAx;
    double *pAz;
    double *Ax;
    double *Ay;
    double *Az;
    double *pAy_gpu;
    double *pAx_gpu;
    double *pAz_gpu;
    double *Energy_gpu;
    cufftDoubleComplex *wfc;
    cufftDoubleComplex *V_gpu;
    cufftDoubleComplex *EV_opt;
    cufftDoubleComplex *wfc_backup;
    cufftDoubleComplex *GK;
    cufftDoubleComplex *GV;
    cufftDoubleComplex *GpAx;
    cufftDoubleComplex *GpAy;
    cufftDoubleComplex *GpAz;
    cufftDoubleComplex *EV;
    cufftDoubleComplex *EK;
    cufftDoubleComplex *EpAy;
    cufftDoubleComplex *EpAx;
    cufftDoubleComplex *EpAz;
    cufftDoubleComplex *EappliedField; 
    cufftDoubleComplex *wfc_gpu;
    cufftDoubleComplex *K_gpu;
    cufftDoubleComplex *par_sum;

    //std::cout << omegaX << '\t' << omegaY << '\n';
    //std::cout << "xDim is: " << xDim << '\t' <<  "yDim is: " << yDim << '\t'
    //          << "zDim is: " << zDim << '\n';

    cufftResult result;
    cufftHandle plan_1d;
    cufftHandle plan_3d;

    dim3 grid = cupar.dim3val("grid");

    std::string buffer;
    double Rxy; //Condensate scaling factor.
    double a0x, a0y, a0z; //Harmonic oscillator length in x and y directions

    unsigned int xD=1,yD=1,zD=1;
    threads = 128;

    // number of blocks in simulation
    unsigned int b = xDim*yDim*zDim/threads;

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
    
    int i, j, k; //Used in for-loops for indexing
    
/*
    double xOffset, yOffset;
    xOffset=0.0;//5.0e-6;
    yOffset=0.0;//5.0e-6;
*/
    
    double mass = 1.4431607e-25; //Rb 87 mass, kg
    par.store("mass",mass);
    double a_s = 4.76e-9;
    par.store("a_s",a_s);

    double sum = 0.0;

    a0x = sqrt(HBAR/(2*mass*omegaX));
    a0y = sqrt(HBAR/(2*mass*omegaY));
    a0z = sqrt(HBAR/(2*mass*omegaZ));
    par.store("a0x",a0x);
    par.store("a0y",a0y);
    par.store("a0z",a0z);

    //std::cout << "a0x and y are: " << a0x << '\t' << a0y << '\n';

    //std::cout << N << '\t' << a_s << '\t' << mass << '\t' << omegaZ << '\n';
    
    Rxy = pow(15,0.2)*pow(N*a_s*sqrt(mass*omegaZ/HBAR),0.2);
    par.store("Rxy",Rxy);
    double bec_length = sqrt( HBAR/(mass*sqrt( omegaX*omegaX * 
                                               ( 1 - omega*omega) ) ));

    //std::cout << "Rxy is: " << Rxy << '\n';
    double xMax = 6*Rxy*a0x; //10*bec_length; //6*Rxy*a0x;
    double yMax = 6*Rxy*a0y; //10*bec_length;
    double zMax = 6*Rxy*a0z; //10*bec_length
    par.store("xMax",xMax);
    par.store("yMax",yMax);
    par.store("zMax",zMax);

    double pxMax, pyMax, pzMax;
    pxMax = (PI/xMax)*(xDim>>1);
    pyMax = (PI/yMax)*(yDim>>1);
    pzMax = (PI/zMax)*(zDim>>1);
    par.store("pyMax",pyMax);
    par.store("pxMax",pxMax);
    par.store("pzMax",pzMax);
    
    double dx = xMax/(xDim>>1);
    double dy = yMax/(yDim>>1);
    double dz = zMax/(zDim>>1);
    par.store("dx",dx);
    par.store("dy",dy);
    par.store("dz",dz);
    
    double dpx, dpy, dpz;
    dpx = PI/(xMax);
    dpy = PI/(yMax);
    dpz = PI/(zMax);
    //std::cout << "yMax is: " << yMax << '\t' << "xMax is: " << xMax << '\n';
    //std::cout << "dpx and dpy are:" << '\n';
    //std::cout << dpx << '\t' << dpy << '\n';
    par.store("dpx",dpx);
    par.store("dpy",dpy);
    par.store("dpz",dpz);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    //double *x,*y,*xp,*yp;
    x = (double *) malloc(sizeof(double) * xDim);
    y = (double *) malloc(sizeof(double) * yDim);
    z = (double *) malloc(sizeof(double) * zDim);
    xp = (double *) malloc(sizeof(double) * xDim);
    yp = (double *) malloc(sizeof(double) * yDim);
    zp = (double *) malloc(sizeof(double) * zDim);

    /*
     * R-space and K-space grids
     */
    //std::cout << "dx and dy are: " << '\n';
    //std::cout << dx << '\t' << dy << '\n';
    // creating x,y,z,xp,yp,zp
    for(i=0; i<xDim/2; ++i){
        x[i] = -xMax + i*dx;        
        x[i + (xDim/2)] = i*dx;
        
        xp[i] = i*dpx;
        xp[i + (xDim/2)] = -pxMax + i*dpx;
        
    }
    for(i=0; i<yDim/2; ++i){
        y[i] = -yMax + i*dy;        
        y[i + (yDim/2)] = i*dy;
        
        yp[i] = i*dpy;
        yp[i + (yDim/2)] = -pyMax + i*dpy;

    }
    for(i=0; i<zDim/2; ++i){
        z[i] = -zMax + i*dz;        
        z[i + (zDim/2)] = i*dz;
        
        zp[i] = i*dpz;
        zp[i + (zDim/2)] = -pzMax + i*dpz;

    }

    par.store("x", x);
    par.store("y", y);
    par.store("z", z);
    par.store("xp", xp);
    par.store("yp", yp);
    par.store("zp", zp);
    

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
    GpAx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    GpAy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    GpAz = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EK = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EV = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EV_opt = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    pAy = (double *) malloc(sizeof(double) * gSize);
    pAx = (double *) malloc(sizeof(double) * gSize);
    pAz = (double *) malloc(sizeof(double) * gSize);
    Ax = (double *) malloc(sizeof(double) * gSize);
    Ay = (double *) malloc(sizeof(double) * gSize);
    Az = (double *) malloc(sizeof(double) * gSize);
    EpAy = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EpAx = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EpAz = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gSize);
    EappliedField = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * 
                                                         gSize);
    
    /* Initialise wfc, EKp, and EVr buffers on GPU */
    cudaMalloc((void**) &Energy_gpu, sizeof(double) * gSize);
    cudaMalloc((void**) &wfc_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &Phi_gpu, sizeof(double) * gSize);
    cudaMalloc((void**) &K_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &V_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &pAy_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &pAx_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &pAz_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void**) &par_sum, sizeof(cufftDoubleComplex) * (gSize/threads));
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //std::cout << "all variables malloc'd" << '\n';

    #ifdef __linux
    int cores = omp_get_num_procs();
    int index;
    par.store("Cores_Total",cores);

    // Assuming dev system specifics (Xeon with HT -> cores detected / 2)
    par.store("Cores_Max",cores/2);
    omp_set_num_threads(cores/2);
    //std::cout << "GAMMAY IS: " << gammaY << '\n';
    //#pragma omp parallel for private(k)
    #endif
    for( i=0; i < xDim; i++ ){
        for( j=0; j < yDim; j++ ){
            for( k=0; k < zDim; k++ ){
                index = i * yDim * zDim + j * zDim + k;
                Phi[index] = fmod(l*atan2(y[j], x[i]),2*PI);
                
                wfc[index].x = exp(-( pow((x[i])/(Rxy*a0x),2) + 
                                      pow((y[j])/(Rxy*a0y),2) +
                                      pow((z[k])/(Rxy*a0z),2)) ) *
                                      cos(Phi[index]);
                wfc[index].y = -exp(-( pow((x[i])/(Rxy*a0x),2) + 
                                       pow((y[j])/(Rxy*a0y),2) +
                                       pow((z[k])/(Rxy*a0z),2)) ) *
                                       sin(Phi[index]);
                
                V[index] = opr.V_fn(par.Vfn)(par, opr, i, j, k);
                K[index] = opr.K_fn(par.Kfn)(par, opr, i, j, k);
    
                GV[index].x = exp( -V[index]*(gdt/(2*HBAR)));
                GK[index].x = exp( -K[index]*(gdt/HBAR));
                GV[index].y = 0.0;
                GK[index].y = 0.0;
    
                // Ax and Ay will be calculated here but are used only for
                // debugging. They may be needed later for magnetic field calc
                Ax[index] = opr.Ax_fn(par.Afn)(par, opr, i, j, k);
                Ay[index] = opr.Ay_fn(par.Afn)(par, opr, i, j, k);
                Az[index] = opr.Az_fn(par.Afn)(par, opr, i, j, k);
                
                pAy[index] = opr.pAy_fn("rotation")(par, opr, i, j, k);
                pAx[index] = opr.pAx_fn("rotation")(par, opr, i, j, k);
                pAz[index] = opr.pAz_fn("rotation")(par, opr, i, j, k);
    
                GpAx[index].x = exp(-pAx[index]*gdt);
                GpAx[index].y = 0;
                GpAy[index].x = exp(-pAy[index]*gdt);
                GpAy[index].y = 0;
                GpAz[index].x = exp(-pAz[index]*gdt);
                GpAz[index].y = 0;
                
                EV[index].x=cos( -V[index]*(dt/(2*HBAR)));
                EV[index].y=sin( -V[index]*(dt/(2*HBAR)));
                EK[index].x=cos( -K[index]*(dt/HBAR));
                EK[index].y=sin( -K[index]*(dt/HBAR));
                
                EpAy[index].x=cos(-omega*omegaX*pAy[index]*dt);
                EpAy[index].y=sin(-omega*omegaX*pAy[index]*dt);
                EpAx[index].x=cos(-omega*omegaX*pAx[index]*dt);
                EpAx[index].y=sin(-omega*omegaX*pAx[index]*dt);
                EpAz[index].x=cos(-omega*omegaX*pAz[index]*dt);
                EpAz[index].y=sin(-omega*omegaX*pAz[index]*dt);
        
                sum+=sqrt(wfc[index].x*wfc[index].x + 
                          wfc[index].y*wfc[index].y);
            }
        }
    }

    std::cout << "writing initial variables to file..." << '\n';
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    //hdfWriteDouble(xDim, V, 0, "V_0"); //HDF COMING SOON!
    //hdfWriteComplex(xDim, wfc, 0, "wfc_0");
    //FileIO::writeOutDouble(buffer, data_dir + "V_opt",V_opt,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "V",V,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "K",K,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "pAy",pAy,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "pAx",pAx,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "pAz",pAz,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "Ax",Ax,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "Ay",Ay,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "Az",Az,gSize,0);
    FileIO::writeOut(buffer, data_dir + "WFC",wfc,gSize,0);
    FileIO::writeOut(buffer, data_dir + "EpAy",EpAy,gSize,0);
    FileIO::writeOut(buffer, data_dir + "EpAx",EpAx,gSize,0);
    FileIO::writeOut(buffer, data_dir + "EpAz",EpAz,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "Phi",Phi,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "r",r,gSize,0);
    FileIO::writeOutDouble(buffer, data_dir + "x",x,xDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "y",y,yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "z",z,zDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "px",xp,xDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "py",yp,yDim,0);
    FileIO::writeOutDouble(buffer, data_dir + "pz",zp,zDim,0);
    FileIO::writeOut(buffer, data_dir + "GK",GK,gSize,0);
    FileIO::writeOut(buffer, data_dir + "GV",GV,gSize,0);
    FileIO::writeOut(buffer, data_dir + "GpAx",GpAx,gSize,0);
    FileIO::writeOut(buffer, data_dir + "GpAy",GpAy,gSize,0);
    FileIO::writeOut(buffer, data_dir + "GpAz",GpAz,gSize,0);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //std::cout << "wrote initial variables" << '\n';

    //free(V); 
    free(K); free(r); free(Phi);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    sum=sqrt(sum*dx*dy*dz);
    //#pragma omp parallel for reduction(+:sum) private(j)
    for (i = 0; i < xDim; i++){
        for (j = 0; j < yDim; j++){
            for (k = 0; k < zDim; k++){
                index = i * yDim * zDim + j * zDim + k;
                wfc[index].x = (wfc[index].x)/(sum);
                wfc[index].y = (wfc[index].y)/(sum);
            }
        }
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    
    //std::cout << "xDim is: " << xDim << '\t' << "yDim is: " << yDim << '\n';
    //std::cout << "plan_2d is: " << plan_2d << '\n';
    result = cufftPlan3d(&plan_3d, xDim, yDim, zDim, CUFFT_Z2Z);
    //std::cout << "found result" << '\n';
    if(result != CUFFT_SUCCESS){
        printf("Result:=%d\n",result);
        printf("Error: Could not execute cufftPlan3d(%s ,%d, %d, %d).\n", 
                "plan_3d",
                (unsigned int)xDim, (unsigned int)yDim, (unsigned int) zDim);
        return -1;
    }

    cufftHandle plan_dim2 = generate_plan_other3d(par, 2);
    cufftHandle plan_dim3 = generate_plan_other3d(par, 3);

    result = cufftPlan1d(&plan_1d, xDim, CUFFT_Z2Z, yDim);
    if(result != CUFFT_SUCCESS){
        printf("Result:=%d\n",result);
        printf("Error: Could not execute cufftPlan1d(%s ,%d ,%d , %d).\n", 
               "plan_1d", (unsigned int)xDim, (unsigned int)yDim,
                          (unsigned int)zDim);
        return -1;
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //std::cout << GV[0].x << '\t' << GK[0].x << '\t' 
    //          << pAy[0] << '\t' << pAx[0] << '\n';

    //std::cout << "storing variables..." << '\n';

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
    par.store("yMax", yMax);
    par.store("winding", l);
    par.store("x", x);
    par.store("y", y);
    par.store("z", z);
    par.store("xp", xp);
    par.store("yp", yp);
    par.store("zp", zp);
    wave.store("Energy", Energy);
    wave.store("r", r);
    opr.store("V", V);
    opr.store("V_opt", V_opt);
    wave.store("Phi", Phi);
    wave.store("Phi_gpu", Phi_gpu);
    opr.store("K", K);
    opr.store("pAy", pAy);
    opr.store("pAx", pAx);
    opr.store("pAz", pAz);
    opr.store("Energy_gpu", Energy_gpu);
    par.store("atoms", N);
    par.store("threads", threads);
    wave.store("wfc", wfc);
    opr.store("V_gpu", V_gpu);
    opr.store("EV_opt", EV_opt);
    wave.store("wfc_backup", wfc_backup);
    opr.store("GK", GK);
    opr.store("GV", GV);
    opr.store("GpAx", GpAx);
    opr.store("GpAy", GpAy);
    opr.store("GpAz", GpAz);
    opr.store("EV", EV);
    opr.store("EK", EK);
    opr.store("EpAy", EpAy);
    opr.store("EpAx", EpAx);
    opr.store("EpAz", EpAz);
    opr.store("EappliedField", EappliedField);
    wave.store("wfc_gpu", wfc_gpu);
    opr.store("K_gpu", K_gpu);
    opr.store("pAy_gpu", pAy_gpu);
    opr.store("pAx_gpu", pAx_gpu);
    opr.store("pAz_gpu", pAz_gpu);
    wave.store("par_sum", par_sum);

    cupar.store("result", result);
    cupar.store("plan_1d", plan_1d);
    cupar.store("plan_3d", plan_3d);
    cupar.store("plan_dim2", plan_dim2);
    cupar.store("plan_dim3", plan_dim3);

    cupar.store("grid", grid);

    std::cout << "variables stored" << '\n';

    return 0;
}

int main(int argc, char **argv){

    Grid par = parseArgs(argc,argv);
    Wave wave;
    Op opr;
    Cuda cupar;

    int device = par.ival("device");
    int dimnum = par.ival("dimnum");
    cudaSetDevice(device);

    std::string buffer;
    time_t start,fin;
    time(&start);
    printf("Start: %s\n", ctime(&start));

    //************************************************************//
    /*
    * Initialise the Params data structure to track params and variables
    */
    //************************************************************//

    // Initialization split between 2d and 3d
    if (dimnum == 2){
        init_2d(opr, cupar, par, wave);
    }
    else if (dimnum == 3){
        init_3d(opr, cupar, par, wave);
    }

    //std::cout << "initialized" << '\n';

    // Re-establishing variables from parsed Grid class
    // Note that 3d variables are set to nullptr's unless needed 
    //      This might need to be fixed later
    std::string data_dir = par.sval("data_dir");
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double *V_opt = opr.dsval("V_opt");
    double *pAy = opr.dsval("pAy");
    double *pAx = opr.dsval("pAx");
    double *pAy_gpu = opr.dsval("pAy_gpu");
    double *pAx_gpu = opr.dsval("pAx_gpu");
    double *pAz_gpu = nullptr;
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    bool read_wfc = par.bval("read_wfc");
    int gsteps = par.ival("gsteps");
    int esteps = par.ival("esteps");
    cufftDoubleComplex *wfc = wave.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *V_gpu = opr.cufftDoubleComplexval("V_gpu");
    cufftDoubleComplex *GK = opr.cufftDoubleComplexval("GK");
    cufftDoubleComplex *GV = opr.cufftDoubleComplexval("GV");
    cufftDoubleComplex *GpAx = opr.cufftDoubleComplexval("GpAx");
    cufftDoubleComplex *GpAy = opr.cufftDoubleComplexval("GpAy");
    cufftDoubleComplex *GpAz = nullptr;
    cufftDoubleComplex *EV = opr.cufftDoubleComplexval("EV");
    cufftDoubleComplex *EK = opr.cufftDoubleComplexval("EK");
    cufftDoubleComplex *EpAy = opr.cufftDoubleComplexval("EpAy");
    cufftDoubleComplex *EpAx = opr.cufftDoubleComplexval("EpAx");
    cufftDoubleComplex *EpAz = nullptr;
    cufftDoubleComplex *wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
    cufftDoubleComplex *K_gpu = opr.cufftDoubleComplexval("K_gpu");
    cufftDoubleComplex *par_sum = wave.cufftDoubleComplexval("par_sum");
    cudaError_t err = cupar.cudaError_tval("err");
    int gsize = xDim * yDim;

    // Special variables for the 3d case
    if (dimnum == 3){
        double dz = par.dval("dz");
        double *z = par.dsval("z");
        double *pAz = opr.dsval("pAz");
        double *pAz_gpu = opr.dsval("pAz_gpu");
        int zDim = par.ival("zDim");
        cufftDoubleComplex *GpAz = opr.cufftDoubleComplexval("GpAz");
        cufftDoubleComplex *EpAz = opr.cufftDoubleComplexval("EpAz");
        gsize = xDim*yDim*zDim;
    }

    std::cout << "variables re-established" << '\n';
    //std::cout << read_wfc << '\n';

    //************************************************************//
    /*
    * Groundstate finder section
    */
    //************************************************************//
    FileIO::writeOutParam(buffer, par, data_dir + "Params.dat");

    // Note: This only works in 2d case
    if(read_wfc){
        printf("Loading wavefunction...");
        wfc=FileIO::readIn("wfc_load","wfci_load",xDim, yDim);
        printf("Wavefunction loaded.\n");
    }

    //std::cout << "gsteps: " << gsteps << '\n';
    
    if(gsteps > 0){
        err=cudaMemcpy(K_gpu, GK, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy K_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(V_gpu, GV, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy V_gpu to device" << '\n';
            exit(1);
        }
        FileIO::writeOut(buffer, data_dir + "GK1",GK,gsize,0);
        FileIO::writeOut(buffer, data_dir + "GV1",GV,gsize,0);
        err=cudaMemcpy(pAy_gpu, GpAy, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy pAy_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(pAx_gpu, GpAx, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy pAx_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy wfc_gpu to device" << '\n';
            exit(1);
        } 
        opr.store("pAx", pAx);
        opr.store("pAy", pAy);
        opr.store("GK", GK);
        opr.store("GV", GV);
        wave.store("wfc", wfc);
        opr.store("K_gpu", K_gpu);
        opr.store("V_gpu", V_gpu);
        wave.store("wfc_gpu", wfc_gpu);
        opr.store("pAy_gpu", pAy_gpu);
        opr.store("pAx_gpu", pAx_gpu);

        // Special cases for 3d
        if (dimnum == 3){

            err=cudaMemcpy(pAz_gpu, GpAz, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            opr.store("pAz_gpu", pAz_gpu);
        
            evolve_3d(wave, opr, par_sum,
                      gsteps, cupar, 0, par, buffer);
        }
        if (dimnum == 2){
            evolve_2d(wave, opr, par_sum,
                      gsteps, cupar, 0, par, buffer);
        }
        wfc = wave.cufftDoubleComplexval("wfc");
        wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
        cudaMemcpy(wfc, wfc_gpu, sizeof(cufftDoubleComplex)*gsize, 
                   cudaMemcpyDeviceToHost);
    }

    std::cout << GV[0].x << '\t' << GK[0].x << '\t' 
              << pAy[0] << '\t' << pAx[0] << '\n';

    //free(GV); free(GK); free(pAy); free(pAx);

    // Re-initializing wfc after evolution
    //wfc = wave.cufftDoubleComplexval("wfc");
    //wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");

    std::cout << "evolution started..." << '\n';
    std::cout << "esteps: " << esteps << '\n';

    //************************************************************//
    /*
    * Evolution
    */
    //************************************************************//
    if(esteps > 0){
        err=cudaMemcpy(pAy_gpu, EpAy, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy pAy_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(pAx_gpu, EpAx, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy pAx_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(K_gpu, EK, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy K_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy V_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*gsize,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy wfc_gpu to device" << '\n';
            exit(1);
        }

        opr.store("pAx", pAx);
        opr.store("pAy", pAy);
        opr.store("EK", EK);
        opr.store("EV", EV);
        wave.store("wfc", wfc);
        opr.store("K_gpu", K_gpu);
        opr.store("V_gpu", V_gpu);
        wave.store("wfc_gpu", wfc_gpu);
        opr.store("pAy_gpu", pAy_gpu);
        opr.store("pAx_gpu", pAx_gpu);

        FileIO::writeOutDouble(buffer, data_dir + "V_opt",V_opt,gsize,0);

        // Special variables / instructions for 3d case
        if (dimnum == 3){
            err=cudaMemcpy(pAz_gpu, EpAz, sizeof(cufftDoubleComplex)*gsize,
                           cudaMemcpyHostToDevice);
            opr.store("pAz_gpu", pAz_gpu);
            evolve_3d(wave, opr, par_sum,
                      esteps, cupar, 1, par, buffer);
        }
        if (dimnum == 2){
            evolve_2d(wave, opr, par_sum,
                      esteps, cupar, 1, par, buffer);
        }
    
        wfc = wave.cufftDoubleComplexval("wfc");
        wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
    }

    std::cout << "done evolving" << '\n';
    free(EV); free(EK); free(EpAy); free(EpAx);
    free(x);free(y);
    cudaFree(wfc_gpu); cudaFree(K_gpu); cudaFree(V_gpu); cudaFree(pAx_gpu); 
    cudaFree(pAy_gpu); cudaFree(par_sum);

    time(&fin);
    printf("Finish: %s\n", ctime(&fin));
    printf("Total time: %ld seconds\n ",(long)fin-start);
    return 0;
}
