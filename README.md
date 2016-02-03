###############################################################################
###############################################################################
 
                                    (
                           (        )\ )               )
                           )\ )    (()/(   (    (     ((
                          (()/(     /(_))  )\   )\   /))
                          /(_))_   (_))   ((_) ((_) ()/(_)
                          ██████╗ ██████╗ ██╗   ██╗███████╗
                         ██╔════╝ ██╔══██╗██║   ██║██╔════╝
                         ██║  ███╗██████╔╝██║   ██║█████╗
                         ██║   ██║██╔═══╝ ██║   ██║██╔══╝
                         ╚██████╔╝██║     ╚██████╔╝███████╗
                          ╚═════╝ ╚═╝      ╚═════╝ ╚══════╝

                    GPUE - GPU Gross-Pitaevskii Equation solver
                             Version 0.abs(exp(i*pi))

                        https://bitbucket.org/loriordan/gpue
                        https://github.com/mlxd/gpue

###############################################################################
## 0. Introduction
###############################################################################

Welcome to GPUE, the [fastest zero temperature BEC routines in the land](http://peterwittek.com/gpe-comparison.html) (the last time I checked).

Runs on CUDA 7.0 (C++11 functionality needed) on both Linux and Mac OS X 
(Nvidia GPU only). I have not tested on Windows. Other requirements are Python 
2.6+ (though PyPy is MUCH faster), Numpy, Scipy, Matplotlib, Mencoder.

To build, first check the predefined paths in the Makefile (CUDA lib/lib64, 
bin, include, etc).

Configuring the simulation parameters is as easy as editing the 
bin/run_params.conf file, and entering them as specified. Each new line will 
represent another simulation, with the maximum number of simultaneous
simulations to be given in run.sh.

To run the simulations:
chmod +x ./run.sh; ./run.sh

If the image generation fails, this can be rectified manually by editing the 
file py/vis.py, and selecting the appropriate operation. Not all Python code 
may currently not behave as expected, due to the switch over to vortex tracking
entirely in C++.

Many routines are also being converted from C host-code to C++11 host-code, 
so if there is any strange behaviour then please let me know. Bug reports are 
always welcome, as well as comments, improvements, cake, free-dinners and cold 
weather.

If you enjoy/use/learn from this code, please give me a citation, as it took a 
long time to develop, and I'd appreciate it greatly! A nice email would suffice 
too :)

###############################################################################
## 1. What does this software do?
###############################################################################

This software is a CUDA-enabled non-linear Schrodinger (Gross-Pitaevskii) 
equation solver. The primary use of this code was for my research on 
rapidly rotating Bose-Einstein condensates. Due to the complexity and 
timescales needed to simulate such system, it was essential to write some 
accelerated code to understand the behaviour of such systems. If you would like 
to know more, check my Google scholar profile 
<https://scholar.google.com/citations?user=s-6jND0AAAAJ&hl=en>
which will update as I finish some papers.

The premise is the following:
We want to simulate how a Bose-Einstein condensate (BEC) behaves in a trap. 
The trap is parabolic (harmonic), and for the lowest energy state of the 
system (ground-state) the BEC will want to sit about the centre. Due to the
interaction between the particles it will occupy more space than a standard 
Schrodinger equation, which has zero interactions. As a result of these 
interactions many interesting things happen.

The main purpose of the code is to investigate the behaviour of quantum 
vortices (think really small tornadoes). Instead of having a continuous 
range of angular momentum values, the condensate can only accept angular 
momentum in quantised predefined units. 

The most interesting fact is that instead of getting bigger and bigger with 
faster rotation (as a tornado would), these vortices only allow themselves 
to enter with a singular unit of angular momentum (think 100x 1 unit vortices 
instead of 1x 100 unit vortex). This gives us a nice well arranged lattice if 
performed correctly. It is this lattice that I have been researching (read as: 
playing with). However, this code can be used in any trapping geometry, 
rotation, etc. that you wish to use. Bear in mind, I have developed this for 
a 2D system only. Extensions to 3D should be straight-forward though.

###############################################################################
## 2. Great! How do I make a BEC?
###############################################################################

Well, first you need to look at the run_params.conf file, and give it some 
necessary parameters to generate your favourite condensate.

As an example
./gpue -x 256 -y 256 -T 1e-3 -t 1e-3 -n 1e5 -g 1e3 -e 1e3 -p 1 -r 0 -w 0.3 -o 10 -d 0 -l 1 -s 1 -i 1.0 -P 0.0 -G 1.0 -L 0 -X 1 -Y 1 -k 0 -O 0.0 -W 1 -U 0 -V 0 -S 0.0

As an example, here are some simulations performed with the code:
https://www.youtube.com/playlist?list=PLiRboSbbz10s6cXxvYLFOn3QbmQpdtQVd

The above parameters will be better explained when the papers have been 
published (which is only fair). Comments will be added for relevant sections 
too to allow for their use. 

###############################################################################
## 3. What do you see as the future of this?
###############################################################################
I would like this tool to be a suite for 1D, 2D and 3D simulations of both 
Schrodinger and non-linear Schrodinger systems. 
