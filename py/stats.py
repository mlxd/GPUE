'''
stats.py - GPUE: Split Operator based GPU solver for Nonlinear 
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
'''
import os
from numpy import genfromtxt
import math as m
#import matplotlib as mpl
import numpy as np
import numpy.matlib
#mpl.use('Agg')
#import multiprocessing as mp
#from multiprocessing import Pool
#from multiprocessing import Process
#from matplotlib.ticker import ScalarFormatter
#import matplotlib.pyplot as plt
import ConfigParser 
import random as r
from decimal import *

#getcontext().prec = 4
c = ConfigParser.ConfigParser()
c.readfp(open(r'Params.dat'))

incr = int(c.getfloat('Params','print_out'))
xDim = int(c.getfloat('Params','xDim'))
yDim = int(c.getfloat('Params','yDim'))

def lsFit(start,end,incr):
	L = np.matrix([
			[0,0,1],
			[1,0,1],
			[0,1,1],
			[1,1,1]
			])
	LSQ = np.linalg.inv(np.transpose(L)*L)*np.transpose(L)
	for i in range(start,end,incr):
		v_arr=genfromtxt('vort_arr_' + str(i),delimiter=',' )
		real=open('wfc_ev_' + str(i)).read().splitlines()
		img=open('wfc_evi_' + str(i)).read().splitlines()
		a_r = np.asanyarray(real,dtype='f8') #64-bit double
		a_i = np.asanyarray(img,dtype='f8') #64-bit double
		a = a_r[:] + 1j*a_i[:]
		wfc = (np.reshape(a,(xDim,yDim)))

		indX = [row[0] for row in v_arr]
		indY = [row[1] for row in v_arr]
		wind = [row[2] for row in v_arr]
		sign = [row[3] for row in v_arr]
		data=[]
		for ii in range(0,len(indX)):
			p=np.matrix([[0],[0],[0],[0]],dtype=np.complex)
			p[0]=(wfc[indX[ii], indY[ii]])
			p[1]=(wfc[indX[ii]+1, indY[ii]])
			p[2]=(wfc[indX[ii], indY[ii]+1])
			p[3]=(wfc[indX[ii]+1, indY[ii]+1])
			rc = LSQ * np.real(p)
			ic = LSQ * np.imag(p)

			A=np.squeeze([row[0:2] for row in [rc,ic]])
			B=-np.squeeze([row[2] for row in [rc,ic]])
			r=np.linalg.lstsq(A,B)[0]
			data.append([indX[ii]+r[0],indY[ii]+r[1],sign[ii]])

#		f = plt.imshow(abs(wfc)**2)
#		plt.jet()
#		plt.gca().invert_yaxis()
#		plt.hold(True)
#		X = [row[0] for row in data]
#		Y = [row[1] for row in data]
#		plt.scatter(Y,X,s=0.2,marker='.',c='red',lw=0)
#		plt.scatter(indY,indX,s=0.2,marker='.',c='yellow',lw=0)
#		plt.savefig("fig.png",dpi=1200)
#		plt.close()
		np.savetxt('vort_lsq_'+str(i)+'.csv',data,delimiter=',')
