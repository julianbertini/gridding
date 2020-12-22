from scipy.io import loadmat
from kb_kernel import calc_kb_kernel
import matplotlib.pyplot as plt
import numpy as np
import grid

if __name__ == "__main__":
    
    mat = loadmat("data/spiralexampledata.mat")
    
    overgrid_factor = 2
    kwidth = 1.5
    klength = 100
    gridsize = 256

    # squeeze mat arrays to remove extra dimensions
    k = np.squeeze(mat['kspacelocations'])
    dcf = np.squeeze(mat['dcf'])
    s = np.squeeze(mat['spiraldata'])

    s_real = s.real.flatten('F')
    s_imag = s.imag.flatten('F')
    sg_real = np.zeros(gridsize**2)
    sg_imag = np.zeros(gridsize**2)
    kx = k.real 
    ky = k.imag 
    kerneltable = calc_kb_kernel(kwidth, overgrid_factor, klength)[0]
    nkernelpoints = kerneltable.size
    nsamples = k.size

    #plt.plot(np.linspace(0, dcf.size, dcf.size), dcf)
    #plt.show()

    #grid.calcdcflut(kx, ky, dcf, kerneltable, gridsize, nsamples, nkernelpoints, kwidth)  

    grid.gridlut(kx, ky, s_real, s_imag, dcf, sg_real, sg_imag, kerneltable, nsamples, gridsize, nkernelpoints, kwidth)
   
    sg = [complex(sg_real[i], sg_imag[i]) for i in range(sg_real.size)]

    sg = np.asarray(sg).reshape(256,256, order='C') 

    im = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(sg)))


    plt.imshow(abs(im))
    plt.show()




