'''
observables.py - GPUE: Split Operator based GPU solver for Nonlinear 
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
import matplotlib as mpl
import numpy as np
import scipy as sp
import numpy.matlib
mpl.use('Agg')
import multiprocessing as mp
from multiprocessing import Pool
from multiprocessing import Process
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import ConfigParser
import random as r
from decimal import *
from scipy.spatial import Delaunay

HBAR = 1.05457148e-34
PI = 3.141592653589793

getcontext().prec = 4
c = ConfigParser.ConfigParser()
c.readfp(open(r'Params.dat'))

xDim = int(c.getfloat('Params','xDim'))
yDim = int(c.getfloat('Params','yDim'))
gndMaxVal = int(c.getfloat('Params','gsteps'))
evMaxVal = int(c.getfloat('Params','esteps'))
incr = int(c.getfloat('Params','print_out'))
#sep = (c.getfloat('Params','dx'))
dx = (c.getfloat('Params','dx'))
dy = (c.getfloat('Params','dy'))
dkx = (c.getfloat('Params','dpx'))
dky = (c.getfloat('Params','dpy'))
dt = (c.getfloat('Params','dt'))
xMax = (c.getfloat('Params','xMax'))
yMax = (c.getfloat('Params','yMax'))
omegaZ = (c.getfloat('Params','omegaZ'))
mass = (c.getfloat('Params','Mass'))
omega = (c.getfloat('Params','omega'))
omegaX = (c.getfloat('Params','omegaX'))

try:
	num_vort = int(c.getfloat('Params','Num_vort'))
except:
	print '!num_vort undefined!'
N = int(c.getfloat('Params','atoms'))

data = numpy.ndarray(shape=(xDim,yDim))

x=np.asarray(open('x_0').read().splitlines(),dtype='f8')
y=np.asarray(open('y_0').read().splitlines(),dtype='f8')
#kx=np.asarray(open('px_0').read().splitlines(),dtype='f8')
#ky=np.asarray(open('py_0').read().splitlines(),dtype='f8')

kx = np.reshape( np.array( [np.linspace( 0, (xDim/2-1)*dkx, xDim/2), np.linspace( (-xDim/2-1)*dkx, -dkx, xDim/2)]), (xDim,1) )
ky = np.reshape( np.array( [np.linspace( 0, (yDim/2-1)*dky, yDim/2), np.linspace( (-yDim/2-1)*dky, -dky, yDim/2)]), (yDim,1) )
kxm, kym = np.meshgrid(kx,ky)
km_mag = np.sqrt( kxm**2 + kym**2 )
k_mag = np.sqrt( kx**2 + ky**2 )
kMax = max(max(k_mag))

hbar = 1.05457e-34 
m = 1.4431607e-25

## Kinetic energy spectrum = kinertrum. Calculates the spectrum for compressible and incompressible kinetic energies.
# @param Psi The wavefunction
# @param dx Increment along x
# @param i The current step number
# @param quOn Boolean to turn on quantum kinetic energy spectrum (includes phase term).
def kinertrum(Psi, dx, i, quOn):

    kMax = np.max(np.max(kx))
    Psi[np.where(Psi==0)] = 1e-100
    n_r = np.abs(Psi)**2
    n_r[np.where(n_r==0)] = 1e-100
    cPsi = np.conj(Psi)
    phi = np.angle(Psi)

    ph1 = np.unwrap(phi, axis=0)
    ph2 = np.unwrap(phi, axis=1)

    vel_ph1_x, vel_ph1_y = np.gradient(ph1,dx,dy)
    vel_ph2_x, vel_ph2_y = np.gradient(ph2,dx,dy)

    v_x = (hbar/m)*vel_ph1_x;
    v_y = (hbar/m)*vel_ph2_y;
    v_x[np.where(v_x==0)] = 1e-100
    v_y[np.where(v_y==0)] = 1e-100

    u_x = np.multiply(np.abs(Psi),v_x)
    u_y = np.multiply(np.abs(Psi),v_y)

    if quOn:
    	u_x = np.multiply(u_x,np.exp(1j*np.angle(Psi)))
    	u_y = np.multiply(u_y,np.exp(1j*np.angle(Psi)))

    u_kx = np.fft.fftn(u_x)
    u_ky = np.fft.fftn(u_y)

    uc_kx = ( kxm**2*u_kx + kxm*kym*u_ky ) / ( km_mag**2 + 1e-100 )
    uc_ky = ( kym*kxm*u_kx + kym**2*u_ky ) / ( km_mag**2 + 1e-100 )

    ui_kx = u_kx - uc_kx
    ui_ky = u_ky - uc_ky

    uc_x = np.fft.ifftn(uc_kx)
    uc_y = np.fft.ifftn(uc_ky)
    ui_x = np.fft.ifftn(ui_kx)
    ui_y = np.fft.ifftn(ui_ky)

    Ec = 0.5*np.abs(np.square(uc_x) + np.square(uc_y))
    Ei = 0.5*np.abs(np.square(ui_x) + np.square(ui_y))

    fig, ax = plt.subplots()
    f = plt.imshow((Ec),cmap=plt.get_cmap('gnuplot2'))
    cbar = fig.colorbar(f)
    plt.gca().invert_yaxis()
    plt.savefig("Ec_" + str(i/incr) + ".png",dpi=200)
    plt.close()
    fig, ax = plt.subplots()
    f = plt.imshow((Ei),cmap=plt.get_cmap('gnuplot2'))
    cbar = fig.colorbar(f)
    plt.gca().invert_yaxis()
    plt.savefig("Ei_" + str(i/incr) + ".png",dpi=200)
    plt.close()
    
    print Ec
    #exit()
    ekc = np.zeros((xDim/2-1,1))
    eki = np.zeros((xDim/2-1,1))
    for i1 in np.arange(0,np.size(k_mag)/2 -2):
        iX = np.array(np.where(np.logical_and( k_mag[i1] >= km_mag, k_mag[i1+1] < km_mag )))
#        Ei_kx = np.sum(np.sum(np.abs(ui_kx[iX]**2*k[iX]))
#        Ei_ky = np.sum(np.sum(np.abs(ui_ky[iX]**2*k[iX]))
        ekc[i1] = (0.5*m*k_mag[i1]) * (np.sum(np.abs(uc_kx[iX]**2 + uc_ky[iX]**2)))/len(iX)
        eki[i1] = (0.5*m*k_mag[i1]) * (np.sum(np.abs(ui_kx[iX]**2 + ui_ky[iX]**2)))/len(iX)
	print i1
    np.savetxt('ekc_' + str(i) + '.csv',ekc,delimiter=',')
    np.savetxt('eki_' + str(i) + '.csv',eki,delimiter=',')
    fig, ax = plt.subplots()
    print eki[0:(xDim/2-2)]
    f = plt.loglog(np.ravel(k_mag[0:(xDim/2 -2)]),eki[0:(xDim/2-2)])
    plt.savefig("eki_" + str(i) + ".png",dpi=200)
    f = plt.loglog(np.ravel(k_mag[0:(xDim/2 -2)]),np.ravel(ekc[0:(xDim/2-2)]))
    plt.savefig("ekc_" + str(i) + ".png",dpi=200)
    plt.close()


def kinertrum_loop(dataName, initValue, finalValue, incr):
	for i in range(initValue,incr*(finalValue/incr),incr):
		if os.path.exists(dataName + '_' + str(i)):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = numpy.asanyarray(real,dtype='f8') #64-bit double
			a_i = numpy.asanyarray(img,dtype='f8') #64-bit double
			a = a_r[:] + 1j*a_i[:]

			kinertrum(np.reshape(a,(xDim,yDim)),dx,i,1)

def dens_struct_fact(dataName, initValue, finalValue,incr):
	n_k=np.zeros(finalValue/incr)
	n_k_t=np.zeros((finalValue/incr,xDim,yDim),dtype=np.complex128)
	for i in range(initValue,incr*(finalValue/incr),incr):
		if os.path.exists(dataName + '_' + str(i)):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = numpy.asanyarray(real,dtype='f8') #64-bit double
			a_i = numpy.asanyarray(img,dtype='f8') #64-bit double
			a = a_r[:] + 1j*a_i[:]
			n = np.abs(a)**2

			kinertrum(np.reshape(a,(xDim,yDim)),dx,i,0)
			sf = np.fft.fftshift(np.fft.fft2(np.reshape(n,(xDim,yDim))))
			n_k_t[i/incr][:][:] = sf[:][:];
			n_k[i/incr]=(abs(np.sum(np.sum(sf))*dkx**2))

			fig, ax = plt.subplots()
			f = plt.imshow(np.log10(abs(sf)),cmap=plt.get_cmap('gnuplot2'))
			cbar = fig.colorbar(f)
			plt.gca().invert_yaxis()
			plt.savefig("struct_" + str(i/incr) + ".png",vmin=0,vmax=12,dpi=200)
			plt.close()
			print i/incr

	np.savetxt('Struct' + '.csv',n_k,delimiter=',')
	plt.plot(range(initValue,finalValue,incr),n_k)
	sp.io.savemat('Struct_t.mat',mdict={'n_k_t',n_k_t})
	plt.savefig("Struct.pdf",dpi=200)
	plt.close()

V = np.array(open('V_0').read().splitlines(),dtype='f8')
V = np.reshape(V,(xDim,yDim))
K = np.array(open('K_0').read().splitlines(),dtype='f8')
K = np.reshape(K,(xDim,yDim))
xPy = np.array(open('xPy_0').read().splitlines(),dtype='f8')
xPy = np.reshape(xPy,(xDim,yDim))
yPx = np.array(open('yPx_0').read().splitlines(),dtype='f8')
yPx = np.reshape(yPx,(xDim,yDim))
g = (0.5*N)*4.0*HBAR*HBAR*PI*(4.67e-9/mass)*np.sqrt(mass*omegaZ/(2.0*PI*HBAR))

def energy_total(dataName, initValue, finalValue, increment):
	E=np.zeros((finalValue,1))
	E_k=np.zeros((finalValue,1))
	E_vi=np.zeros((finalValue,1))
	E_l=np.zeros((finalValue,1))
	for i in range(initValue,incr*(finalValue/incr),incr):
		if os.path.exists(dataName + '_' + str(i)):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = np.array(real,dtype='f8') #64-bit double
			a_i = np.array(img,dtype='f8') #64-bit double
			wfcr = np.reshape(a_r[:] + 1j*a_i[:],(xDim,yDim))
			wfcp = np.array(np.fft.fft2(wfcr))
			wfcr_c = np.conj(wfcr)
			
			E1 = np.fft.ifft2(K*wfcp)
			E2 = (V + 0.5*g*np.abs(wfcr)**2)*wfcr
			E3 = -(omega*omegaX)*(np.fft.ifft(xPy*np.fft.fft(wfcr,axis=0),axis=0) - np.fft.ifft(yPx*np.fft.fft(wfcr,axis=1),axis=1)  )
			
			E_k[i/incr] = np.trapz(np.trapz(wfcr_c*E1))*dx*dy
			E_vi[i/incr] = np.trapz(np.trapz(wfcr_c*E2))*dx*dy
			E_l[i/incr] = np.trapz(np.trapz(wfcr_c*E3))*dx*dy
			E[i/incr] = E_k[i/incr] + E_vi[i/incr] + E_l[i/incr]
			print (i/float(evMaxVal))
	np.savetxt('E_'+ str(i) + '.csv',E,delimiter=',')
	np.savetxt('E_k_'+ str(i) + '.csv',E_k,delimiter=',')
	np.savetxt('E_vi_'+ str(i) + '.csv',E_vi,delimiter=',')
	np.savetxt('E_l_'+ str(i) + '.csv',E_l,delimiter=',')
	t = np.array(range(initValue,finalValue,incr))/dt
	plt.plot(t,E,'r-',t,E_k,'g-',t,E_vi,'b-',t,E_l,'y-')
	plt.savefig("EnergyVst.pdf",dpi=200)
	plt.close()

def energy_kinetic(dataName, initValue, finalValue, increment):
	px1 = np.fft.fftshift(px)
	py1 = np.fft.fftshift(py)
	dk=[]
	dk2[:] = (px1[:]**2 + py1[:]**2)
	Lz = np.zeros( (finalValue/incr))
	for i in range(initValue,incr*(finalValue/incr),incr):
		if os.path.exists(dataName + '_' + str(i)):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = numpy.asanyarray(real,dtype='f8') #64-bit double
			a_i = numpy.asanyarray(img,dtype='f8') #64-bit double
			a = a_r[:] + 1j*a_i[:]
			wfcp = np.fft.fft2(np.reshape(a,(xDim,yDim)))
			conjwfcp = np.conj(wfcp)
			E_k = np.zeros(len(px1))
			for ii in range(0,len(px1)):
				E_k[ii] = np.sum( np.sum( np.multiply(wfcp,conjwfcp) )  )*dk2[ii]

		np.savetxt('E_k_' + str(i) + '.csv',E_k,delimiter=',')
		print i

def energy_potential(dataName, initValue, finalValue, increment):
	print 'energy'

def ang_mom(dataName, initValue, finalValue, incr, ev_type, imgdpi):
	xm, ym = np.meshgrid(x,y)
	pxm, pym = np.meshgrid(px,py)
	dx2=dx**2
	Lz = np.zeros( (finalValue/incr))
	for i in range(initValue,incr*(finalValue/incr),incr):
		if os.path.exists(dataName + '_' + str(i)):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = numpy.asanyarray(real,dtype='f8') #64-bit double
			a_i = numpy.asanyarray(img,dtype='f8') #64-bit double
			a = a_r[:] + 1j*a_i[:]
			wfc = np.reshape(a,(xDim,yDim))
			conjwfc = np.conj(wfc)

			wfc_ypx = np.multiply(ym,np.fft.ifft(np.multiply(pxm,np.fft.fft(wfc,axis=1)),axis=1))
			wfc_xpy = np.multiply(xm,np.fft.ifft(np.multiply(pym,np.fft.fft(wfc,axis=0)),axis=0))
			result = np.sum( np.sum( np.multiply(conjwfc,wfc_xpy - wfc_ypx) ) )*dx2
		else:
			print "Skipped " + dataName + "_"+ str(i)
			result = np.nan

		print i, incr
		Lz[(i/incr)] = np.real(result)
	type=""
	if ev_type == 0:
		type = "gnd"
	else:
		type = "ev"
	np.savetxt('Lz.csv',Lz,delimiter=',')

	plt.plot(Lz)
	plt.savefig("Lz_"+type+".pdf",dpi=imgdpi)
	plt.axis('off')
	plt.savefig("Lz_"+type+"_axis0.pdf",bbox_inches='tight',dpi=imgdpi)
	plt.close()

def expec_val_monopole(dataName, initValue, finalValue, incr):
	x=np.asarray(open('x_0').read().splitlines(),dtype='f8')
	y=np.asarray(open('y_0').read().splitlines(),dtype='f8')
#	px=open('px_0')
#	py=open('py_0')
	xm, ym = np.meshgrid(x, y)
	result = []
	for i in range(initValue,finalValue,incr):
		if not os.path.exists(dataName):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = numpy.asanyarray(real,dtype='f8') #64-bit double
			a_i = numpy.asanyarray(img,dtype='f8') #64-bit double
			a = a_r[:] + 1j*a_i[:]
			wfc = np.reshape(a,(xDim,yDim))
			conjwfc = np.conj(wfc)

			d1 = np.multiply( np.square(xm) + np.square(ym), wfc )
			d2 = np.multiply( conjwfc, d1)
			result.append( np.real( np.sum( np.sum( d2  ) ) )*dx*dx )
		print str(100*float(i)/finalValue) + '%'
	np.savetxt('monopole.csv',result,delimiter=',')
	plt.plot(range(initValue,finalValue,incr),result)
	plt.savefig("Monopole.png",dpi=200)
	plt.close()

def expec_val_quadrupole(dataName, initValue, finalValue, incr):
	x=np.asarray(open('x_0').read().splitlines(),dtype='f8')
	y=np.asarray(open('y_0').read().splitlines(),dtype='f8')
#	px=open('px_0')
#	py=open('py_0')
	xm, ym = np.meshgrid(x, y)
	result = []
	for i in range(initValue,finalValue,incr):
		if not os.path.exists(dataName):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = numpy.asanyarray(real,dtype='f8') #64-bit double
			a_i = numpy.asanyarray(img,dtype='f8') #64-bit double
			a = a_r[:] + 1j*a_i[:]
			wfc = np.reshape(a,(xDim,yDim))
			conjwfc = np.conj(wfc)

			d1 = np.multiply( np.square(xm) - np.square(ym), wfc )
			d2 = np.multiply( conjwfc, d1)
			result.append( np.real( np.sum( np.sum( d2  ) ) )*dx*dx )
		print str(100*float(i)/finalValue) + '%'
	np.savetxt('quadrupole.csv',result,delimiter=',')
	plt.plot(range(initValue,finalValue,incr),result)
	plt.savefig("Quadrupole.png",dpi=200)
	plt.close()

def expec_val_(quant_name, quantity, dataName, initValue, finalValue, incr):
	x=np.asarray(open('x_0').read().splitlines(),dtype='f8')
	y=np.asarray(open('y_0').read().splitlines(),dtype='f8')
#	px=open('px_0')
#	py=open('py_0')
	xm, ym = np.meshgrid(x, y)
	result = []
	for i in range(initValue,finalValue,incr):
		if not os.path.exists(dataName):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = numpy.asanyarray(real,dtype='f8') #64-bit double
			a_i = numpy.asanyarray(img,dtype='f8') #64-bit double
			a = a_r[:] + 1j*a_i[:]
			wfc = np.reshape(a,(xDim,yDim))
			conjwfc = np.conj(wfc)

			d1 = np.multiply( quantity, wfc )
			d2 = np.multiply( conjwfc, d1)
			result.append( np.real( np.sum( np.sum( d2  ) ) )*dx*dx )
		print str(100*float(i)/finalValue) + '%'
	np.savetxt(quant_name + '.csv',result,delimiter=',')
	plt.plot(range(initValue,finalValue,incr),result)
	plt.savefig(quant_name + ".pdf",dpi=200)
	plt.close()

if __name__ == '__main__':
	kinertrum_loop('wfc_ev', 0, evMaxVal, incr)
	exit()
	energy_total('wfc_ev',0,evMaxVal,incr)
	dens_struct_fact('wfc_ev', 0, evMaxVal, 500)
	
	energy_kinetic('wfc_ev', 0, evMaxVal, 200)
#	ang_mom('wfc_0_ramp', 0, gndMaxVal, incr, 0, 200)
	ang_mom('wfc_ev', 0, evMaxVal, incr, 1, 200)
	expec_val_monopole('wfc_ev',0,evMaxVal,incr)
	expec_val_quadrupole('wfc_ev',0,evMaxVal,incr)
