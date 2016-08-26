# Simple script to plot variables. Will grow with time.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse

parser = argparse.ArgumentParser(description='reading strings for plotting')
parser.add_argument('strings', metavar='ID', nargs='+', help='string to plot')

args = parser.parse_args()
print(args.strings)

#class for all variables with initial definitions
class params:

    # Defining static / default values
    xDim = 512
    yDim = 512
    
    # data_dir is assumed to be in the previous directory
    data_dir = "data"
    
    # Defaulting to first element after imaginary time evolution
    start = 0
    end = 1
    incr = 1

    # item to work with
    item = "wfc"

# Function to plot specific variable
def plot_var(xDim, yDim, data_dir, pltval):
    data = "../" + data_dir + "/" + pltval
    lines = np.loadtxt(data)
    val = np.reshape(lines, (xDim,yDim))
    plt.imshow(val, extent=(1,xDim,1,yDim), interpolation='nearest',
                   cmap = cm.jet)
    plt.colorbar()
    plt.show()
    
# Function to plot wfc with pltvar as a variable to modify the type of plot
def plot_wfc(xDim, yDim, data_dir, pltval, start, end, incr):
    for i in range(start,end,incr):
        data_real = "../" + data_dir + "/wfc_0_const_%s" % i
        data_im = "../" + data_dir + "/wfc_0_consti_%s" % i
    
        lines_real = np.loadtxt(data_real)
        lines_im = np.loadtxt(data_im)
        wfc_real = np.reshape(lines_real, (xDim,yDim));
        wfc_im = np.reshape(lines_im, (xDim,yDim));
    
        wfc = abs(wfc_real + 1j * wfc_im)
        wfc = wfc * wfc
    
        #wfc_k = np.fft.fft2(wfc) 
        #wfc_k_plot = np.abs(np.fft.fftshift(wfc_k))
        #wfc_k_plot = wfc_k_plot**2
        
        plt.imshow(wfc, extent=(1,xDim,1,yDim), interpolation='nearest',
                   cmap = cm.jet)
        plt.colorbar()
        plt.show()
        #fig = plt.figure()
        #fig.savefig('wfc.png')

# Function to plot complex vals with pltvar as the variable
def plot_complex(xDim, yDim, data_dir, pltval, start, end, incr):

    data_real = "../" + data_dir + "/" + pltval + "_0"
    data_im = "../" + data_dir + "/" + pltval + "i_0"
    
    lines_real = np.loadtxt(data_real)
    lines_im = np.loadtxt(data_im)
    wfc_real = np.reshape(lines_real, (xDim,yDim));
    wfc_im = np.reshape(lines_im, (xDim,yDim));
    
    wfc = abs(wfc_real + 1j * wfc_im)
    wfc = wfc * wfc
        
    plt.imshow(wfc, extent=(1,xDim,1,yDim), interpolation='nearest',
               cmap = cm.jet)
    plt.colorbar()
    plt.show()
    #fig = plt.figure()
    #fig.savefig('wfc.png')

# Function to plot wfc with pltvar as a variable to modify the type of plot
def plot_wfc_k(xDim, yDim, data_dir, pltval, start, end, incr):
    for i in range(start,end,incr):
        data_real = "../" + data_dir + "/wfc_0_const_%s" % i
        data_im = "../" + data_dir + "/wfc_0_consti_%s" % i

        lines_real = np.loadtxt(data_real)
        lines_im = np.loadtxt(data_im)
        wfc_real = np.reshape(lines_real, (xDim,yDim));
        wfc_im = np.reshape(lines_im, (xDim,yDim));

        wfc = (wfc_real + 1j * wfc_im)

        wfc_k = np.fft.fft2(wfc)
        wfc_k_plot = np.abs(np.fft.fftshift(wfc_k))
        wfc_k_plot = wfc_k_plot**2

        plt.imshow(wfc_k_plot, extent=(1,xDim,1,yDim), interpolation='nearest',
                   cmap = cm.jet)
        plt.colorbar()
        plt.show()
        #fig = plt.figure()
        #fig.savefig('wfc.png')

# Function to plot wfc with pltvar as a variable to modify the type of plot
def plot_wfc_phase(xDim, yDim, data_dir, pltval, start, end, incr):
    for i in range(start,end,incr):
        data_real = "../" + data_dir + "/wfc_0_const_%s" % i
        data_im = "../" + data_dir + "/wfc_0_consti_%s" % i

        lines_real = np.loadtxt(data_real)
        lines_im = np.loadtxt(data_im)
        wfc_real = np.reshape(lines_real, (xDim,yDim));
        wfc_im = np.reshape(lines_im, (xDim,yDim));

        wfc = (wfc_real + 1j * wfc_im)

        wfc = np.angle(wfc)
        plt.imshow(wfc, extent=(1,xDim,1,yDim), interpolation='nearest',
                   cmap = cm.jet)
        plt.colorbar()
        plt.show()
        #fig = plt.figure()
        #fig.savefig('wfc.png')


# Function to parse arguments for plotting
# Note: We assume that the parameters come in sets
def parse_args(string_list):
    i = 0
    par = params()
    print(string_list[i])
    while i < len(string_list):
        # -d for "data_dir"
        if (string_list[i] == "d"):
            par.data_dir = string_list[i+1]
            i += 2
        # -i for "item" -- The thing to plot
        elif (string_list[i] == "i"):
            par.item = string_list[i+1]
            i += 2
        # -r for "range"
        elif (string_list[i] == "r"):
            par.first = int(string_list[i+1])
            par.end = int(string_list[i+2])
            par.incr = int(string_list[i+3])
            i+= 4
        # -g for "grid"
        elif (string_list[i] == "g"):
            par.xDim = int(string_list[i+1])
            par.yDim = int(string_list[i+2])
            i += 3
    return par

def plot(par):
    if (par.item == "wfc"):
        plot_wfc(par.xDim, par.yDim, par.data_dir, par.item, 
                 par.start, par.end, par.incr)
    elif (par.item == "wfc_k"):
        plot_wfc_k(par.xDim, par.yDim, par.data_dir, par.item,
                   par.start, par.end, par.incr)
    elif (par.item == "wfc_phase"):
        plot_wfc_phase(par.xDim, par.yDim, par.data_dir, par.item,
                       par.start, par.end, par.incr)
    elif (par.item == "GK" or par.item == "GV"):
        plot_complex(par.xDim, par.yDim, par.data_dir, par.item,
                     par.start, par.end, par.incr)
    else:
        plot_var(par.xDim, par.yDim, par.data_dir, par.item)

par = parse_args(args.strings)
plot(par)
