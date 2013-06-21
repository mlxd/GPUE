#
# visual_ev.py - GPUE: Split Operator based GPU solver for Nonlinear 
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

import numpy, scipy, mayavi
from mayavi.mlab import *
from scipy.io import *
import numpy.matlib
xDim=256
yDim=256
data = numpy.ndarray(shape=(xDim,yDim))
s = "./wfc"
figure(size=(xDim,yDim))
for i in range(0,50000,500):
	real=open(s + '_' + str(i)).read().splitlines()
	img=open(s + 'i_' + str(i)).read().splitlines()
	a_r = numpy.asanyarray(real,dtype='f8') #128-bit complex
	a_i = numpy.asanyarray(img,dtype='f8') #128-bit complex
	a = a_r[:] + 1j*a_i[:]
	b = numpy.reshape(a,(xDim,yDim))
	imshow(abs(b)**2)
	view(0,0)
	savefig("wfc_"+str(i)+".png")#,size=(800,600))
#	close(gcf())	
del a, a_r, a_i
#contour3d(b, contours=4, transparent=True)
#imshow(abs(b)**2)
#data_tpot = scipy.io.loadmat('/home/mlxd/workspace/Dev/Tpot.mat')
#oct_a = data_tpot['Pot']
#contour3d(oct_a, contours=4, transparent=True)
#data_wfc = scipy.io.loadmat('/home/mlxd/workspace/Dev/WFC_0.mat')
#oct_b = data_wfc['wfabs']
#contour3d(oct_b, contours=4, transparent=True)
