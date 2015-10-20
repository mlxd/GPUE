'''
vort.py - GPUE: Split Operator based GPU solver for Nonlinear 
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
###############################################################################
import os
from numpy import genfromtxt
import math as m
import numpy as np
import copy as cp
import ConfigParser 

###############################################################################
c = ConfigParser.ConfigParser()
c.readfp(open(r'Params.dat'))

xDim = int(c.getfloat('Params','xDim'))
yDim = int(c.getfloat('Params','yDim'))
gndMaxVal = int(c.getfloat('Params','gsteps'))
evMaxVal = int(c.getfloat('Params','esteps'))
incr = int(c.getfloat('Params','print_out'))
dx = (c.getfloat('Params','dx'))
dt = (c.getfloat('Params','dt'))
xMax = (c.getfloat('Params','xMax'))
yMax = (c.getfloat('Params','yMax'))

###############################################################################
class Vortex: #Tracks indivisual vortices over time.
###############################################################################
###############################################################################
	def __init__(self,uid,x,y,isOn,sign=1):
###############################################################################
		self.uid = uid
		self.x = x
		self.y = y
		self.sign = sign
		self.isOn = isOn
		self.next = None

###############################################################################
	def update_uid(self,uid):
###############################################################################
		self.uid = uid

###############################################################################
	def update_on(self,isOn): #Vortex is trackable
###############################################################################
		self.isOn = isOn

###############################################################################
	def update_next(self,next): #Get next vortex
###############################################################################
		self.next = next

###############################################################################
	def dist(self,vtx): #Distance between self and vtx
###############################################################################
		r = m.sqrt((self.x - vtx.x)**2 + (self.y - vtx.y)**2)
		return r
	
###############################################################################
class VtxList: #Linked-list for tracking vortices
###############################################################################
###############################################################################
	def __init__(self):
###############################################################################
		self.head = None
		self.tail = None
		self.length = 0

###############################################################################
	def element(self,pos): #Get vtx at position pos
###############################################################################
		pos_l = 0
		if pos < self.length:
			vtx = self.head
			while pos_l < pos:
				pos_l = pos_l +1
				vtx = vtx.next
		else:
			print "Out of bounds"
			exit(-1)
		return vtx

###############################################################################
	def vtx_uid(self,uid): #Get vtx with identifier uid
###############################################################################
		vtx = self.head
		pos = 0
		while vtx.uid != uid:
			vtx = vtx.next
			pos = pos +1
		return [vtx,pos]
		
###############################################################################
	def max_uid(self): #Return position and value of largest uid
###############################################################################
		val = 0
		vtx = self.head
		val = vtx.uid
		pos = 0
		#while pos < self.length:
		while True:
			vtx = vtx.next
			if(vtx == None):
				break
			if vtx.uid > val:
				val = vtx.uid
			pos = pos +1
		return [val,pos]
		
###############################################################################
	def add(self,Vtx,index=None): #Add a vtx at index, otherwise end
###############################################################################
		if self.length == 0:
			self.head = Vtx
			self.tail = Vtx
			self.length = 1
		elif index == None:
			self.tail.next = Vtx
			self.tail = Vtx
			self.length = self.length +1
		else:
			Vtx.next = self.element(index)
			self.element(index-1).next = Vtx
			self.length = self.length + 1	
	
###############################################################################
	def as_np(self): #Return numpy array with format x,y,sign,uid,isOn
###############################################################################
		dtype = [('x',float),('y',float),('sign',int),('uid',int),('isOn',int)]
		data =[]# np.array([],dtype=dtype)
		i = 0
		vtx = self.head
		while vtx != None:
			data.append([vtx.x, vtx.y, vtx.sign, vtx.uid, vtx.isOn])
			vtx = vtx.next
			i = i+1
		return (data)

###############################################################################
	def write_out(self,time,data): #Write out CSV file as  x,y,sign,uid,isOn
###############################################################################
		np.savetxt('vort_ord_'+str(time)+'.csv',data,fmt='%10.5f,%10.5f,%i,%i,%i',delimiter=',')

###############################################################################
	def idx_min_dist(self,vortex, isSelf=False): #Closest vtx to self
###############################################################################
		counter = 0
		ret_idx = counter
		vtx = self.head
		if vtx != None:
			r = vtx.dist(vortex)
			while vtx.next != None:
				vtx = vtx.next
				counter = counter +1
				if r > vtx.dist(vortex):
					r = vtx.dist(vortex)
					ret_idx = counter
		return (ret_idx,r)

###############################################################################
	def remove(self,pos): #Remove vortices outside articificial boundary
###############################################################################
		if self.length > 1 and pos > 1:
			current = self.element(pos-1).next
			self.element(pos - 1).next = current.next
			current.next = None
			self.length = self.length - 1
			return current
		elif pos == 0:
			current = self.head
			self.head = self.head.next
			self.length = self.length - 1
			return current
		else:
			self.head = None
			self.length = 0
			return None

###############################################################################
	def swap_uid(self,uid_i,uid_f): #Swap uid between vtx
###############################################################################
		vtx_pos = self.vtx_uid(uid_i)
		self.remove(pos_i)
		self.add(vtx,index=pos_f)

###############################################################################
	def vort_decrease(self,positions,vorts_p): #Turn off vortex timeline
###############################################################################
		max_uid = vorts_p.max_uid()
		for i4 in positions:
			vtx = cp.copy(i4)
			vtx.update_on(False)
			vtx.update_next(None)
			self.add(vtx)
		
###############################################################################
	def vort_increase(self,positions,vorts_p): #Add new vtx to tracking
###############################################################################
		counter = 1
		max_uid = vorts_p.max_uid()
		for i4 in positions:
			self.element(i4).update_uid(max_uid[0] + counter)
			counter = counter+1
		
###############################################################################
def do_the_thing(start,fin,incr): #Performs the tracking
###############################################################################
	#v_arr_p=genfromtxt('vort_lsq_' + str(0) + '.csv',delimiter=',')
	v_arr_p=genfromtxt('vort_arr_' + str(1000),delimiter=',')
	for i in range( start + incr + 1000, fin + 1, incr): #loop over samples in time
	#	print v_arr_p[:,2]
		vorts_p = VtxList()
		vorts_c = VtxList()
		#v_arr_c=genfromtxt('vort_lsq_' + str(i) + '.csv',delimiter=',' )
		v_arr_c=genfromtxt('vort_arr_' + str(i), delimiter=',')
		if i==2000:
			v_arr_p_coords = np.array([a for a in v_arr_p[:,[1,3]]])
			v_arr_c_coords = np.array([a for a in v_arr_c[:,[1,3]]])
			v_arr_p_sign = np.array([a for a in v_arr_p[:,4]])
			v_arr_c_sign = np.array([a for a in v_arr_c[:,4]])
		else:
			v_arr_p_coords = np.array([a for a in v_arr_p[:,[0,1]]])
			v_arr_p_sign = np.array([a for a in v_arr_p[:,2]])
			
		v_arr_c_coords = np.array([a for a in v_arr_c[:,[1,3]]])
		v_arr_c_sign = np.array([a for a in v_arr_c[:,4]])
		for i1 in range(0,v_arr_p_coords.size/2): #loop over coordinates for a given time
			vtx_p = Vortex(i1,v_arr_p_coords[i1][0],v_arr_p_coords[i1][1],True,sign=v_arr_p_sign[i1])#,v_arr_p[i1][2])
			vorts_p.add(vtx_p)
		
		for i2 in range(0,v_arr_c_coords.size/2):
			vtx_c = Vortex(-1-i2,v_arr_c_coords[i2][0],v_arr_c_coords[i2][1],True,sign=v_arr_c_sign[i2])#,v_arr_p[i1][0])
			vorts_c.add(vtx_c)

		for i3 in range(0,vorts_p.length):
			index_r = vorts_c.idx_min_dist(vorts_p.element(i3))
			
			v0c = vorts_c.element(index_r[0]).sign
			v0p = vorts_p.element(i3).sign
			v1c = vorts_c.element(index_r[0]).uid
			if (index_r[1] < 7) and (vorts_c.element(index_r[0]).sign == vorts_p.element(i3).sign) and (vorts_c.element(index_r[0]).uid < 0):
			#if (index_r[1] < 2) and (vorts_c.element(index_r[0]).sign > 0) and (vorts_c.element(index_r[0]).uid < 0):
				vorts_c.element(index_r[0]).update_uid(vorts_p.element(i3).uid)
				vorts_c.element(index_r[0]).update_on(True)

		#You will never remember why this works
		uid_c = [[a for a in b][3] for b in vorts_c.as_np()]
		uid_p = [[a for a in b][3] for b in vorts_p.as_np()]
		#Check the difference between current and previous vtx data
		dpc = set(uid_p).difference(set(uid_c))
		dcp = set(uid_c).difference(set(uid_p))
		vtx_pos_p=[]
		vtx_pos_c=[]
		for i5 in dpc:
			vtx_pos_p = np.append(vtx_pos_p,vorts_p.vtx_uid(i5)[0])
		for i6 in dcp:	
			vtx_pos_c = np.append(vtx_pos_c,vorts_c.vtx_uid(i6)[1])
		if len(dpc or dcp) >= 1:
			vorts_c.vort_decrease(vtx_pos_p,vorts_p)
			vorts_c.vort_increase(vtx_pos_c,vorts_p)

		vorts_c_update=sorted(vorts_c.as_np(),key=lambda vtx: vtx[3])
		vorts_c.write_out(i,np.asarray(vorts_c_update))
		print "[" + str(i) +"]", "Length of previous=" + str(len(v_arr_p_coords)), "Length of current=" + str(len(vorts_c_update))
		v_arr_p=genfromtxt('vort_ord_' + str(i) + '.csv',delimiter=',' )

###############################################################################
###############################################################################
do_the_thing(0,evMaxVal,incr)
