'''
vis.py - GPUE: Split Operator based GPU solver for Nonlinear 
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
CPUs = os.environ['SLURM_JOB_CPUS_PER_NODE']
print "Number of cores: " + str(CPUs)
from numpy import genfromtxt
import math as m
import matplotlib as mpl
import matplotlib.tri as tri
import numpy as np
import scipy as sp
from scipy.spatial import Voronoi, voronoi_plot_2d
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
import stats
import hist3d
import mpld3
from mpld3 import plugins

getcontext().prec = 4
c = ConfigParser.ConfigParser()
getcontext().prec = 4
c = ConfigParser.ConfigParser()
c.readfp(open(r'Params.dat'))

xDim = int(c.getfloat('Params','xDim'))
yDim = int(c.getfloat('Params','yDim'))
gndMaxVal = int(c.getfloat('Params','gsteps'))
evMaxVal = int(c.getfloat('Params','esteps'))
incr = int(c.getfloat('Params','print_out'))
sep = (c.getfloat('Params','dx'))
dx = (c.getfloat('Params','dx'))
dt = (c.getfloat('Params','dt'))
xMax = (c.getfloat('Params','xMax'))
yMax = (c.getfloat('Params','yMax'))
num_vort = 0#int(c.getfloat('Params','Num_vort'))

data = numpy.ndarray(shape=(xDim,yDim))

def delaunay(dataName,dataType,value):
	v_arr=genfromtxt(dataName + str(value) + dataType,delimiter=',' )
	data = np.array([[row[0],row[1]] for row in v_arr])
	dln = sp.spatial.Delaunay(data)
	plt.triplot(data[:,0],data[:,1],dln.simplices.copy(),linewidth=0.5,color='b',marker='.')
	plt.xlim(300,700);plt.ylim(300,700);
	plt.savefig('delaun_' + str(value) + '.png',dpi=200)
	print 'Saved Delaunay @ t=' + str(value)

def voronoi(dataName,dataType,value):
	v_arr=genfromtxt(dataName + str(value) + dataType,delimiter=',' )
	data = [[row[0],row[1]] for row in v_arr]
	vor = Voronoi(data)
	voronoi_plot_2d(vor)
	plt.xlim(300,700);plt.ylim(300,700);
	plt.savefig('voronoi_' + str(value) + '.png',dpi=200)
	print 'Saved Voronoi @ t=' + str(value)

def laplacian(density,name,imgdpi):
	gx,gy = np.gradient(density)
	g2x,gxgy = np.gradient(gx)
	gygx,g2y = np.gradient(gy)
	fig, ax = plt.subplots()
	#f = plt.quiver(gx,gy)
	f = plt.imshow((g2x**2 + g2y**2),cmap=plt.get_cmap('spectral'))
	cbar = fig.colorbar(f)
	plt.savefig(name + "_laplacian.png",dpi=imgdpi)
	plt.close()
	f = plt.imshow((gxgy - gygx),cmap=plt.get_cmap('spectral'))
	cbar = fig.colorbar(f)
	plt.savefig(name + "_dxdy.png",dpi=imgdpi)
	plt.close()

def struct_fact(density,name,imgdpi):
	fig, ax = plt.subplots()
	#f = plt.quiver(gx,gy)
	f = plt.imshow((np.abs(np.fft.fftshift(np.fft.fft2(density)))),cmap=plt.get_cmap('prism'))
	cbar = fig.colorbar(f)
	cbar.set_clim(1e6,1e11)
	plt.jet()
	plt.savefig(name + "_struct_log10.png",dpi=imgdpi)
	plt.close()

def opPot(dataName,imgdpi):
	data = open(dataName).read().splitlines()
	a = numpy.asanyarray(data,dtype='f8')
	b = np.reshape(a,(xDim,yDim))
	fig, ax = plt.subplots()
	f = plt.imshow((b))
	plt.gca().invert_yaxis()
	cbar = fig.colorbar(f)
	plt.jet()
	plt.savefig(dataName + ".png",dpi=imgdpi)
	plt.close()

def hist_gen(name,value,num_bins):
	v_arr=genfromtxt('vort_arr_' + str(value),delimiter=',' )
	H=[]
	count=0

	for i1 in range(0,v_arr.size/2):
		for i2 in range(i1,v_arr.size/2):
			H.append(m.sqrt( abs(v_arr[i1][0]*sep - v_arr[i2][0]*sep)**2  +  abs(v_arr[i1][1]*sep - v_arr[i2][1]*sep)**2 ))
			count = count + 1
	plt.title('Vortex lattice @ t=' + str(value*dt))
	plt.ticklabel_format(style='scientific')
	plt.ticklabel_format(style='scientific',axis='x', scilimits=(0,0))
	h = plt.hist(H, bins=num_bins)
	plt.savefig(name + "_" + str(value) + ".pdf")
	plt.close()

def image_gen(dataName, initValue, finalValue, increment,imgdpi):
	for i in range(initValue,finalValue,increment):
		if not os.path.exists(dataName+"r_"+str(i)+"_abspsi2.png"):
			real=open(dataName + '_' + str(i)).read().splitlines()
			img=open(dataName + 'i_' + str(i)).read().splitlines()
			a_r = numpy.asanyarray(real,dtype='f8') #64-bit double
			a_i = numpy.asanyarray(img,dtype='f8') #64-bit double
			a = a_r[:] + 1j*a_i[:]
			b = np.reshape(a,(xDim,yDim))
			f = plt.imshow(abs(b)**2)
			plt.jet()
			plt.gca().invert_yaxis()
			plt.savefig(dataName+"r_"+str(i)+"_abspsi2.png",dpi=imgdpi)
			plt.close()
			g = plt.imshow(np.angle(b))
			plt.gca().invert_yaxis()
			plt.savefig(dataName+"r_"+str(i)+"_phi.png",dpi=imgdpi)
			plt.close()
			f = plt.imshow(abs(np.fft.fftshift(np.fft.fft2(b)))**2)
			plt.gca().invert_yaxis()
			plt.jet()
			plt.savefig(dataName+"p_"+str(i)+"_abspsi2.png",dpi=imgdpi)
			plt.close()
			g = plt.imshow(np.angle(np.fft.fftshift(np.fft.fft2(b))))
			plt.gca().invert_yaxis()
			plt.savefig(dataName+"p_"+str(i)+"_phi.png",dpi=imgdpi)
			plt.close()
			print "Saved figure: " + str(i) + ".png"
			plt.close()
		else:
			print "File(s) " + str(i) +".png already exist."

def image_gen_single(dataName, value, imgdpi,opmode):
	real=open(dataName + '_' + str(0)).read().splitlines()
	img=open(dataName + 'i_' + str(0)).read().splitlines()
	a1_r = numpy.asanyarray(real,dtype='f8') #128-bit complex
	a1_i = numpy.asanyarray(img,dtype='f8') #128-bit complex
	a1 = a1_r[:] + 1j*a1_i[:]
	b1 = np.reshape(a1,(xDim,yDim))

	if not os.path.exists(dataName+"r_"+str(value)+"_abspsi2.png"):
		real=open(dataName + '_' + str(value)).read().splitlines()
		img=open(dataName + 'i_' + str(value)).read().splitlines()
		a_r = numpy.asanyarray(real,dtype='f8') #128-bit complex
		a_i = numpy.asanyarray(img,dtype='f8') #128-bit complex
		a = a_r[:] + 1j*a_i[:]
		b = np.reshape(a,(xDim,yDim))
		m_val=np.max(np.abs(b)**2)
		#scaleAxis(b,dataName,"_abspsi2",value,imgdpi)
		if opmode & 0b100000 > 0:
#			fig, ax = plt.subplots()
#			#plt.rc('text',usetex=True)
#			#plt.rc('font',family='serif')
#			f = plt.imshow((abs(b)**2 - abs(b1)**2),cmap='gnuplot2',vmin=-6,vmax=6)
#			plt.title(r'$\left(\rho( r,t ) - \rho( r,t_0 )\right),t=$' + str(value*dt))
#			cbar = fig.colorbar(f)
#			plt.gca().set_xlabel('x '+ str((dx)))
#			plt.gca().set_ylabel('x '+ str(dx))
#			plt.gca().invert_yaxis()
#			plt.savefig(dataName+"r_"+str(value)+"_diffabspsi2.png",dpi=imgdpi)
#			plt.close()
#			#plt.rc('text',usetex=True)
#			#plt.rc('font',family='serif')

			fig, ax = plt.subplots()
			f = plt.imshow((abs(b)**2),cmap='gnuplot2',vmin=0,vmax=1e7)
			plt.title('rho(r) @ t=' + str(value*dt))
		#	plt.title(r'$\\rho \left( r,t \right),\,t=$' + str(value*dt))
			
			#plugins.connect(fig, plugins.MousePosition(fontsize=14))
			
			cbar = fig.colorbar(f)
			plt.gca().set_xlabel('x '+ str((dx)))
			plt.gca().set_ylabel('x '+ str(dx))
			plt.gca().invert_yaxis()
			plt.savefig(dataName+"r_"+str(value)+"_abspsi2.png",dpi=imgdpi)
			plt.axis('off')
			plt.savefig(dataName+"r_"+str(value)+"_abspsi2_axis0.pdf",bbox_inches='tight',dpi=imgdpi)
			plt.close()

		if opmode & 0b010000 > 0:
			fig, ax = plt.subplots()
			g = plt.imshow(np.angle(b))
			cbar = fig.colorbar(g)
			plt.gca().invert_yaxis()
			plt.title('theta(r) @ t=' + str(value*dt))
			plt.savefig(dataName+"r_"+str(value)+"_phi.png",dpi=imgdpi)
			plt.close()

		if opmode & 0b001000 > 0:
			fig, ax = plt.subplots()
			f = plt.imshow(abs(np.fft.fftshift(np.fft.fft2(b)))**2)
			cbar = fig.colorbar(f)
			plt.gca().invert_yaxis()
			plt.jet()
			plt.title('rho(p) @ t=' + str(value*dt))
			plt.savefig(dataName+"p_"+str(value)+"_abspsi2.png",dpi=imgdpi)
			plt.close()

		if opmode & 0b000100 > 0:
			fig, ax = plt.subplots()
			g = plt.imshow(np.angle(np.fft.fftshift(np.fft.fft2(b))))
			cbar = fig.colorbar(g)
			plt.gca().invert_yaxis()
			plt.title('theta(p) @ t=' + str(value*dt))
			plt.savefig(dataName+"p_"+str(value)+"_phi.png",dpi=imgdpi)
			plt.close()

		if opmode & 0b000010 > 0:
			struct_fact(abs(b)**2,dataName+"_" + str(value),imgdpi)

		if opmode & 0b000001 > 0:
			laplacian(abs(b)**2,dataName+"_" + str(value),imgdpi)

		print "Saved figure: " + str(value) + ".png"
		plt.close()
	else:
		print "File(s) " + str(value) +".png already exist."

def vort_traj(name,imgdpi):
	evMaxVal_l = evMaxVal
	H=genfromtxt('vort_arr_0',delimiter=',' )
	count=0
	for i1 in range(incr,evMaxVal_l,incr):
		try:
			v_arr=genfromtxt('vort_lsq_' + str(i1) + '.csv',delimiter=',' )
			H=np.column_stack((H,v_arr))
		except:
			evMaxVal_l = i1
			break
	X=np.zeros((evMaxVal_l/incr),dtype=np.float64)
	Y=np.zeros((evMaxVal_l/incr),dtype=np.float64)
	H=np.reshape(H,([num_vort,2,evMaxVal_l/incr]),order='F')
	for i1 in range(0, num_vort):
		for i2 in range(0,evMaxVal_l/incr):
			X[i2]=(H[i1,0,i2]*dx) - xMax
			Y[i2]=(H[i1,1,i2]*dx) - yMax
		h = plt.plot(X,Y,color=(r.random(),r.random(),r.random(),0.85),linewidth=0.1)
	plt.axis('equal')
	plt.title('Vort(x,y) from t=0 to t='+str(evMaxVal_l*dt)+" s")

	plt.axis((-xMax/2.0, xMax/2.0, -yMax/2.0, yMax/2.0))
	plt.ticklabel_format(style='scientific')
	plt.ticklabel_format(style='scientific',axis='x', scilimits=(0,0))
	plt.ticklabel_format(style='scientific',axis='y', scilimits=(0,0))
	plt.savefig(name +".pdf")
	plt.close()
	print "Trajectories plotted."

def scaleAxis(data,dataName,label,value,imgdpi):
	fig, ax = plt.subplots()
	ax.xaxis.set_major_locator(ScaledLocator(dx=dx))
	ax.xaxis.set_major_formatter(ScaledLocator(dx=dx))
	f = plt.imshow(abs(data)**2)
	cbar = fig.colorbar(f)
	plt.gca().invert_yaxis()
	plt.jet()
	plt.savefig(dataName+"r_"+str(value)+"_"+label +".png",dpi=imgdpi)
	plt.close()

def overlap(dataName, initValue, finalValue, increment):
	real=open(dataName + '_' + str(0)).read().splitlines()
	img=open(dataName + 'i_' + str(0)).read().splitlines()
	a_r = numpy.asanyarray(real,dtype='f8') #128-bit complex
	a_i = numpy.asanyarray(img,dtype='f8') #128-bit complex
	wfc0 = a_r[:] + 1j*a_i[:]
	for i in range(initValue,finalValue,increment):
		real=open(dataName + '_' + str(value)).read().splitlines()
		img=open(dataName + 'i_' + str(value)).read().splitlines()
		a_r = numpy.asanyarray(real,dtype='f8') #128-bit complex
		a_i = numpy.asanyarray(img,dtype='f8') #128-bit complex
		a = a_r[:] + 1j*a_i[:]
		b = np.dot(wfc0,a)
		print i, np.sum(b)

if __name__ == '__main__':
	try:
		delaunay('vort_arr_',0)
		stats.lsFit(0,evMaxVal,incr)
		hist3d.plot_hist_pcolor(0,evMaxVal,incr,'b')
		vort_traj('traj_plot',200)
	except:
		print "Unhandled error occurred. Blame Lee."
	opPot('V_opt_0',200)
	opPot('V_0',200)
	opPot('K_0',200)
	gndImgList=[]
	evImgList=[]
	for i in range(0,gndMaxVal,incr):
		gndImgList.append(i)
	for i in range(0,evMaxVal,incr):
		evImgList.append(i)
	gnd_proc = []
	ev_proc = []
	while gndImgList:
		i=gndImgList.pop()
		gnd_proc.append(Process(target=image_gen_single,args=("wfc_0_ramp",i,200,0b110000)))
		gnd_proc.append(Process(target=image_gen_single,args=("wfc_0_const",i,200,0b110000)))
	while evImgList:
		i=evImgList.pop()
		ev_proc.append(Process(target=image_gen_single,args=("wfc_ev",i,200,0b101000)))
		#ev_proc.append(Process(target=mpld3.show,))
		ev_proc.append(Process(target=delaunay,args=("vort_lsq_",'.csv',i)))
		ev_proc.append(Process(target=voronoi,args=("vort_lsq_",'.csv',i)))
		ev_proc.append(Process(target=hist_gen,args=("hist_ev",i,128)))
	proc = gnd_proc + ev_proc
	while proc:
		#if (mp.cpu_count()/2) > len(mp.active_children()):
		if int(CPUs) > len(mp.active_children()):
			print len(mp.active_children())
			try:
				p=proc.pop()
				p.start()
			except:
				print "Failed to execute ", p
