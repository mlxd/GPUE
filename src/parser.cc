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

// note: read variables into appendData directly.

Grid parseArgs(int argc, char** argv){

    // Creates initial grid class for the parameters
    Grid par;
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
                par.store("xDim",xDim);
                break;
            }
            case 'y':
            {
                int yDim = atoi(optarg);
                printf("Argument for y is given as %d\n",yDim);
                par.store("yDim",yDim);
                break;
            }
            case 'z':
            {
                int zDim = atoi(optarg);
                printf("Argument for z is given as %d\n",zDim);
                par.store("zDim",zDim);
                break;
            }
            case 'w':
            {
                double omega = atof(optarg);
                printf("Argument for OmegaRotate is given as %E\n",omega);
                par.store("omega",omega);
                break;
            }
            case 'G':
            {
                double gammaY = atof(optarg);
                printf("Argument for gamma is given as %E\n",gammaY);
                par.store("gammaY",gammaY);
                break;
            }
            case 'g':
            {
                double gsteps = atof(optarg);
                printf("Argument for Groundsteps is given as %d\n",gsteps);
                par.store("gsteps",gsteps);
                break;
            }
            case 'e':
            {
                double esteps = atof(optarg);
                printf("Argument for EvSteps is given as %d\n",esteps);
                par.store("esteps",esteps);
                break;
            }
            case 'T':
            {
                double gdt = atof(optarg);
                printf("Argument for groundstate Timestep is given as %E\n",
                       gdt);
                par.store("gdt",gdt);
                break;
            }
            case 't':
            {
                double dt = atof(optarg);
                printf("Argument for Timestep is given as %E\n",dt);
                par.store("dt",dt);
                break;
            }
            case 'd':
            {
                int device = atoi(optarg);
                printf("Argument for device is given as %d\n",device);
                par.store("device",device);
                break;
            }
            case 'n':
            {
                double atoms = atof(optarg);
                printf("Argument for atoms is given as %d\n",atoms);
                par.store("atoms",atoms);
                break;
            }
            case 'r':
            {
                int read_wfc  = atoi(optarg);
                printf("Argument for ReadIn is given as %d\n",read_wfc);
                par.store("read_wfc",read_wfc);
                break;
            }
            case 'p':
            {
                int print = atoi(optarg);
                printf("Argument for Printout is given as %d\n",print);
                par.store("print_out",print);
                break;
            }
            case 'L':
            {
                double l = atof(optarg);
                printf("Vortex winding is given as : %E\n",l);
                par.store("winding",l);
                break;
            }
            case 'l':
            {
                int ang_mom = atoi(optarg);
                printf("Angular Momentum mode engaged: %d\n",ang_mom);
                par.store("corotating",ang_mom);
                break;
            }
            case 's':
            {
                int gpe = atoi(optarg);
                printf("Non-linear mode engaged: %d\n",gpe);
                par.store("gpe",gpe);
                break;
            }
            case 'o':
            {
                double omegaZ = atof(optarg);
                printf("Argument for OmegaZ is given as %E\n",omegaZ);
                par.store("omegaZ",omegaZ);
                break;
            }
            case 'i':
            {
                double interaction = atof(optarg);
                printf("Argument for interaction scaling is %E\n",interaction);
                par.store("int_scaling",interaction);
                break;
            }
            case 'P':
            {
                double laser_power = atof(optarg);
                printf("Argument for laser power is %E\n",laser_power);
                par.store("laser_power",laser_power);
                break;
            }
            case 'X':
            {
                double omegaX = atof(optarg);
                printf("Argument for omegaX is %E\n",omegaX);
                par.store("omegaX",omegaX);
                break;
            }
            case 'Y':
            {
                double omegaY = atof(optarg);
                printf("Argument for omegaY is %E\n",omegaY);
                par.store("omegaY",omegaY);
                break;
            }
            case 'O':
            {
                double angle_sweep = atof(optarg);
                printf("Argument for angle_sweep is %E\n",angle_sweep);
                par.store("angle_sweep",angle_sweep);
                break;
            }
            case 'k':
            {
                int kick_it = atoi(optarg);
                printf("Argument for kick_it is %i\n",kick_it);
                par.store("kick_it",kick_it);
                break;
            }
            case 'W':
            {
                int write_it = atoi(optarg);
                printf("Argument for write_it is %i\n",write_it);
                par.store("write_it",write_it);
                break;
            }
            case 'U':
            {
                double x0_shift = atof(optarg);
                printf("Argument for x0_shift is %lf\n",x0_shift);
                par.store("x0_shift",x0_shift);
                break;
            }
            case 'V':
            {
                double y0_shift = atof(optarg);
                printf("Argument for y0_shift is %lf\n",y0_shift);
                par.store("y0_shift",y0_shift);
                break;
            }
            case 'S':
            {
                double sepMinEpsilon = atof(optarg);
                printf("Argument for sepMinEpsilon is %lf\n",sepMinEpsilon);
                par.store("sepMinEpsilon",sepMinEpsilon);
                break;
            }
            case 'a':
            {
                int graph = atoi(optarg);
                printf("Argument for graph is %d\n",graph);
                par.store("graph",graph);
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
                return par;
            default:
                abort ();
            }
        }
    }

    return par;
}
