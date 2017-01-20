#-------------gen_data.py------------------------------------------------------##
# Purpose: This file will take the data from GPUE and turn it into a bvox file
#          for visualization with blender
#
#------------------------------------------------------------------------------#

import numpy as np

xDim = yDim = zDim = 128

# Function to plot wfc with pltvar as a variable to modify the type of plot
def wfc_density(xDim, yDim, zDim, data_dir, pltval, i):
    print(i)
    data_real = "../" + data_dir + "/wfc_0_const_%s" % i
    data_im = "../" + data_dir + "/wfc_0_consti_%s" % i
    if (pltval == "wfc_ev"):
        data_real = "../" + data_dir + "/wfc_ev_%s" % i
        data_im = "../" + data_dir + "/wfc_evi_%s" % i
    lines_real = np.loadtxt(data_real)
    lines_im = np.loadtxt(data_im)
    print(len(lines_real))
    wfc_real = np.reshape(lines_real, (xDim,yDim,zDim));
    wfc_im = np.reshape(lines_im, (xDim,yDim, zDim));
    wfc = abs(wfc_real + 1j * wfc_im)
    wfc = wfc * wfc
    wfc = np.reshape(wfc,(xDim*yDim*zDim))
    maximum = max(wfc)
    wfc /= maximum
    wfc = np.reshape(wfc,(xDim,yDim,zDim))
    #print(wfc)
    return wfc

def wfc_phase(xDim, yDim, zDim, data_dir, pltval, i):
    print(i)
    data_real = "../" + data_dir + "/wfc_0_const_%s" % i
    data_im = "../" + data_dir + "/wfc_0_consti_%s" % i
    if (pltval == "wfc_ev"):
        data_real = "../" + data_dir + "/wfc_ev_%s" % i
        data_im = "../" + data_dir + "/wfc_evi_%s" % i
    lines_real = np.loadtxt(data_real)
    lines_im = np.loadtxt(data_im)
    wfc_real = np.reshape(lines_real, (xDim,yDim, zDim));
    wfc_im = np.reshape(lines_im, (xDim,yDim, zDim));
    wfc = (wfc_real + 1j * wfc_im)
    wfc = np.angle(wfc)
    wfc = np.reshape(wfc,(xDim*yDim*zDim))
    maximum = max(wfc)
    minimum = min(wfc)
    for i in range(0,len(wfc)):
        wfc[i] = (wfc[i] - minimum) / (maximum - minimum)
    wfc = np.reshape(wfc,(xDim,yDim,zDim))
    return wfc

def var(xDim, yDim, zDim, data_dir, pltval):
    data = "../" + data_dir + "/" + pltval
    lines = np.loadtxt(data)
    maximum = max(lines)
    minimum = min(lines)
    for i in range(0,len(lines)):
        lines[i] = (lines[i] - minimum) / (maximum - minimum)
    val = np.reshape(lines, (xDim,yDim,zDim));
    return val

def proj_var2d(xdim, yDim, zDim, data_dir, pltval):
    filename = "../" + data_dir + "/" + "val"
    file = open(filename,"w")
    data = "../" + data_dir + "/" + pltval
    lines = np.loadtxt(data)
    var_data = np.reshape(lines, (xDim, yDim, zDim))
    for k in range(0,xDim):
        for j in range(0,yDim):
            file.write(str(var_data[k][j][zDim / 2])+'\n')
    file.close


def proj_2d(xDim, yDim, zDim, data_dir, pltval, i):
    filename = "../" + data_dir + "/wfc_0"
    print(i)
    data_real = "../" + data_dir + "/wfc_0_const_%s" % i
    data_im = "../" + data_dir + "/wfc_0_consti_%s" % i
    if (pltval == "wfc_ev"):
        data_real = "../" + data_dir + "/wfc_ev_%s" % i
        data_im = "../" + data_dir + "/wfc_evi_%s" % i
    lines_real = np.loadtxt(data_real)
    lines_im = np.loadtxt(data_im)
    print(len(lines_real))
    wfc_real = np.reshape(lines_real, (xDim,yDim,zDim));
    wfc_im = np.reshape(lines_im, (xDim,yDim, zDim));
    wfc = abs(wfc_real + 1j * wfc_im)
    wfc = wfc * wfc
    file = open(filename,'w')
    for k in range(0,xDim):
        for j in range(0,yDim):
            file.write(str(wfc[j][k][zDim/2]) + '\n')
    file.close()

def proj_k2d(xDim, yDim, zDim, data_dir, pltval, i):
    filename = "../" + data_dir + "/wfc_1"
    print(i)
    data_real = "../" + data_dir + "/wfc_0_const_%s" % i
    data_im = "../" + data_dir + "/wfc_0_consti_%s" % i
    lines_real = np.loadtxt(data_real)
    lines_im = np.loadtxt(data_im)
    print(len(lines_real))
    wfc_real = np.reshape(lines_real, (xDim,yDim,zDim));
    wfc_im = np.reshape(lines_im, (xDim,yDim, zDim));
    wfc = np.fft.fftshift(np.fft.fftn(wfc_real + 1j * wfc_im))
    wfc = abs(wfc) * abs(wfc)
    file = open(filename,'w')
    for k in range(0,xDim):
        for j in range(0,yDim):
            file.write(str(wfc[j][k][zDim/2]) + '\n')
    file.close()

#proj_var2d(xDim, yDim, zDim, "data", "Ax_0")
#proj_2d(xDim, yDim, zDim,"data","wfc",50000)
#proj_k2d(xDim, yDim, zDim,"data","wfc",50000)

#item = wfc_density(xDim, yDim, zDim,"data","wfc",50000)
#item = wfc_phase(xDim, yDim, zDim,"data","wfc",90)
item = var(xDim, yDim, zDim,"data","gauge_3d")

nx, ny, nz, nframes = xDim, yDim, zDim,1
header = np.array([nx,ny,nz,nframes])

binfile = open('test1.bvox','wb')
header.astype('<i4').tofile(binfile)
item.astype('<f4').tofile(binfile)
