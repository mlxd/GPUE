/*
* parser.cc - GPUE: Split Operator based GPU solver for Nonlinear 
Schrodinger Equation, Copyright (C) 2011-2016, James Ryan Schloss
<jrs.schloss@gmail.com>, Tadhg Morgan, Neil Crowley. All rights reserved.

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
#include "../include/parser.h"
#include "../include/ds.h"

// note: read variables into appendData directly.
/*

void parseArgs(int argc, char** argv, Array params){

    int opt;
    while ((opt = getopt (argc, argv, 
           "d:x:y:w:G:g:e:T:t:n:p:r:o:L:l:s:i:P:X:Y:O:k:W:U:V:S:a:")) != -1)
    {
        switch (opt)
        {
            case 'x':
            {
                int xDim = atoi(optarg);
                printf("Argument for x is given as %d\n",xDim);
                appendData(&params,"xDim",(double)xDim);
                break;
            }
            case 'y':
            {
                int yDim = atoi(optarg);
                printf("Argument for y is given as %d\n",yDim);
                appendData(&params,"yDim",(double)yDim);
                break;
            }
            case 'w':
            {
                double omega = atof(optarg);
                printf("Argument for OmegaRotate is given as %E\n",omega);
                appendData(&params,"omega",omega);
                break;
            }
            case 'G':
            {
                double gammaY = atof(optarg);
                printf("Argument for gamma is given as %E\n",gammaY);
                appendData(&params,"gammaY",gammaY);
                break;
            }
            case 'g':
            {
                double gsteps = atof(optarg);
                printf("Argument for Groundsteps is given as %ld\n",gsteps);
                appendData(&params,"gsteps",gsteps);
                break;
            }
            case 'e':
            {
                double esteps = atof(optarg);
                printf("Argument for EvSteps is given as %ld\n",esteps);
                appendData(&params,"esteps",esteps);
                break;
            }
            case 'T':
            {
                double gdt = atof(optarg);
                printf("Argument for groundstate Timestep is given as %E\n",
                       gdt);
                appendData(&params,"gdt",gdt);
                break;
            }
            case 't':
            {
                double dt = atof(optarg);
                printf("Argument for Timestep is given as %E\n",dt);
                appendData(&params,"dt",dt);
                break;
            }
            case 'd':
            {
                int device = atoi(optarg);
                printf("Argument for device is given as %d\n",device);
                appendData(&params,"device",device);
                break;
            }
            case 'n':
            {
                double atoms = atof(optarg);
                printf("Argument for atoms is given as %ld\n",atoms);
                appendData(&params,"atoms",atoms);
                break;
            }
            case 'r':
            {
                int read_wfc  = atoi(optarg);
                printf("Argument for ReadIn is given as %d\n",read_wfc);
                appendData(&params,"read_wfc",(double)read_wfc);
                break;
            }
            case 'p':
            {
                int print = atoi(optarg);
                printf("Argument for Printout is given as %d\n",print);
                appendData(&params,"print_out",(double)print);
                break;
            }
            case 'L':
            {
                double l = atof(optarg);
                printf("Vortex winding is given as : %E\n",l);
                appendData(&params,"winding",l);
                break;
            }
            case 'l':
            {
                int ang_mom = atoi(optarg);
                printf("Angular Momentum mode engaged: %d\n",ang_mom);
                appendData(&params,"corotating",(double)ang_mom);
                break;
            }
            case 's':
            {
                int gpe = atoi(optarg);
                printf("Non-linear mode engaged: %d\n",gpe);
                appendData(&params,"gpe",gpe);
                break;
            }
            case 'o':
            {
                double omegaZ = atof(optarg);
                printf("Argument for OmegaZ is given as %E\n",omegaZ);
                appendData(&params,"omegaZ",omegaZ);
                break;
            }
            case 'i':
            {
                double interaction = atof(optarg);
                printf("Argument for interaction scaling is %E\n",interaction);
                appendData(&params,"int_scaling",interaction);
                break;
            }
            case 'P':
            {
                double laser_power = atof(optarg);
                printf("Argument for laser power is %E\n",laser_power);
                appendData(&params,"laser_power",laser_power);
                break;
            }
            case 'X':
            {
                double omegaX = atof(optarg);
                printf("Argument for omegaX is %E\n",omegaX);
                appendData(&params,"omegaX",omegaX);
                break;
            }
            case 'Y':
            {
                double omegaY = atof(optarg);
                printf("Argument for omegaY is %E\n",omegaY);
                appendData(&params,"omegaY",omegaY);
                break;
            }
            case 'O':
            {
                double angle_sweep = atof(optarg);
                printf("Argument for angle_sweep is %E\n",angle_sweep);
                appendData(&params,"angle_sweep",angle_sweep);
                break;
            }
            case 'k':
            {
                int kick_it = atoi(optarg);
                printf("Argument for kick_it is %i\n",kick_it);
                appendData(&params,"kick_it",kick_it);
                break;
            }
            case 'W':
            {
                int write_it = atoi(optarg);
                printf("Argument for write_it is %i\n",write_it);
                appendData(&params,"write_it",write_it);
                break;
            }
            case 'U':
            {
                double x0_shift = atof(optarg);
                printf("Argument for x0_shift is %lf\n",x0_shift);
                appendData(&params,"x0_shift",x0_shift);
                break;
            }
            case 'V':
            {
                double y0_shift = atof(optarg);
                printf("Argument for y0_shift is %lf\n",y0_shift);
                appendData(&params,"y0_shift",y0_shift);
                break;
            }
            case 'S':
            {
                double sepMinEpsilon = atof(optarg);
                printf("Argument for sepMinEpsilon is %lf\n",sepMinEpsilon);
                appendData(&params,"sepMinEpsilon",sepMinEpsilon);
                break;
            }
            case 'a':
            {
                int graph = atoi(optarg);
                printf("Argument for graph is %d\n",graph);
                appendData(&params,"graph",graph);
                break;
            }
            case '?':
            {
                if (optopt == 'c') {
                    fprintf (stderr, 
                             "Option -%c requires an argument.\n", optopt);
                } 
                else if (isprint (optopt)) {
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                } 
                else {
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",optopt);
                }
                return -1;
            default:
                abort ();
            }
        }
    }
}
*/
