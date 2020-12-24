from scipy.io import loadmat
from kb_kernel import calc_kb_kernel
import matplotlib.pyplot as plt
import numpy as np
import grid

def read_klocations(filename):
    """ Reads a .txt containing k-space locations of data samples
    """
    kx_locs = []
    ky_locs = []
    with open(filename, 'r') as f:
        for line in f:
            kx_ky_tokens = line.split(',')
            kx_tokens = kx_ky_tokens[0].split('=')
            kx_locs.append(float(kx_tokens[1].strip()))
            ky_tokens = kx_ky_tokens[1].split('=')
            ky_locs.append(float(ky_tokens[1].strip()))

    return kx_locs, ky_locs

def avg_kpoints(data):
    """ Average every 4 data samples, since MR sequence 
        records 4 ADC points for each k-space location.

        Data should be trimmed before using this method.
    """
    avg_data = np.zeros((int(data.shape[0]/4),data.shape[1],data.shape[2]), dtype=complex)
    for i in range(0,data.shape[0],4):
        j = int(i/4)
        avg_data[j,:,:] = np.average([data[i,:,:],data[i+1,:,:],data[i+2,:,:],data[i+3,:,:]]) 

    return avg_data 

def trim_kpoints(data):
    """ Trim ends of k-space, since they contain ramp-up
        and ramp-down points.

        These points do not correspond to valuable signal.

        Start idx: 480 (120 * 4)
        End idx: 31920 (32000 - 4 * 20) 
    """
    return data[480:31920,:,:]

if __name__ == "__main__":
    
    #mat = loadmat("data/spiralexampledata.mat")
    ball_phantom_mat = loadmat("./data/ball_phantom_data/ball_phantom.mat")

    kx, ky = read_klocations("./data/ball_phantom_data/rosetteKs.txt")
    
    overgrid_factor = 2.37
    kwidth = 1.5
    klength = 500
    gridsize = 256

    # squeeze mat arrays to remove extra dimensions
    #k = np.squeeze(mat['kspacelocations'])
    #dcf = np.squeeze(mat['dcf'])
    dcf = np.ones_like(kx)
    #s = np.squeeze(mat['spiraldata'])
    s = np.squeeze(ball_phantom_mat['data'])
    s = trim_kpoints(s)
    s = avg_kpoints(s)

    # Need to be in col major order (Fortran) to match with the way C code is written
    # This is not intuitive, since we would expect C major order.
    s_real = s.real.flatten('F')
    s_imag = s.imag.flatten('F')
    sg_real = np.zeros(gridsize**2)
    sg_imag = np.zeros(gridsize**2)
    #kx = k.real 
    #ky = k.imag 
    kerneltable = calc_kb_kernel(kwidth, overgrid_factor, klength)[0]
    nkernelpoints = kerneltable.size
    nsamples = len(kx)

    grid.calcdcflut(kx, ky, dcf, kerneltable, gridsize, nsamples, nkernelpoints, kwidth)  
    #plt.plot(np.linspace(0, dcf.size, dcf.size), dcf)
    #plt.show()

    grid.gridlut(kx, ky, s_real, s_imag, dcf, sg_real, sg_imag, kerneltable, nsamples, gridsize, nkernelpoints, kwidth)
   
    sg = [complex(sg_real[i], sg_imag[i]) for i in range(sg_real.size)]

    sg = np.asarray(sg).reshape(256,256, order='F') 

    im = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(sg)))

    plt.imshow(abs(im))
    plt.show()
