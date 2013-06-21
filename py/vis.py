#
# vis.py - GPUE: Split Operator based GPU solver for Nonlinear 
# Schrodinger Equation, Copyright (C) 2012, Lee J. O'Riordan, Tadhg 
# Morgan, Neil Crowley. 

# This library is free software; you can redistribute it and/or modify 
# it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation; either version 2.1 of the 
# License, or (at your option) any later version. This library is 
# distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
# License for more details. You should have received a copy of the GNU 
# Lesser General Public License along with this library; if not, write 
# to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
# Boston, MA 02111-1307 USA 
#
import matplotlib as mpl
import scipy
import numpy as np
from scipy.io import *
import numpy.matlib
mpl.use('Agg')
import matplotlib.pyplot as plt
xDim=2048/8
yDim=2048/8
data = numpy.ndarray(shape=(xDim,yDim))
s = "./wfc_0"
for i in range(0,100000,1000):
	real=open(s + '_' + str(i)).read().splitlines()
	img=open(s + 'i_' + str(i)).read().splitlines()
	a_r = numpy.asanyarray(real,dtype='f8') #128-bit complex
	a_i = numpy.asanyarray(img,dtype='f8') #128-bit complex
	a = a_r[:] + 1j*a_i[:]
	b = np.reshape(a,(xDim,yDim))
	f = plt.imshow(abs(b)**2)
	plt.jet()
	plt.savefig("wfc_"+str(i)+".png",dpi=200)
	plt.close()
	g = plt.imshow(numpy.angle(b))
	plt.savefig("phi_"+str(i)+".png",dpi=200)
	plt.close()
	h = plt.contour(abs(b)**2,80,linewidths=0.5)
	plt.savefig("cntr_"+str(i)+".png",dpi=200)
	print "Saved figure: " + str(i) + ".png"
	plt.close()
s = "./wfc"
for i in range(0,10000,200):
	real=open(s + '_' + str(i)).read().splitlines()
	img=open(s + 'i_' + str(i)).read().splitlines()
	a_r = numpy.asanyarray(real,dtype='f8') #128-bit complex
	a_i = numpy.asanyarray(img,dtype='f8') #128-bit complex
	a = a_r[:] + 1j*a_i[:]
	b = np.reshape(a,(xDim,yDim))
	f = plt.imshow(abs(b)**2)
	plt.jet()
	plt.savefig("wfc_ev_"+str(i)+".png",dpi=200)
	plt.close()
	g = plt.imshow(numpy.angle(b))
	plt.savefig("phi_ev_"+str(i)+".png",dpi=200)
	plt.close()
	h = plt.contour(abs(b)**2,80,linewidths=0.5)
	plt.savefig("cntr_ev_"+str(i)+".png",dpi=200)
	print "Saved figure: " + str(i) + ".png"
	plt.close()
del a, a_r, a_i
