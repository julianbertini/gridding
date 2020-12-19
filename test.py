from scipy.io import loadmat
from kb_kernel import calc_kb_kernel
import matplotlib.pyplot as plt
import numpy as np
import grid

if __name__ == "__main__":
    
    mat = loadmat("./spiralexampledata.mat")

    overgrid_factor = 2
    kwidth = 1.5
    klength = 100
    gridsize = 256

    # squeeze mat arrays to remove extra dimensions
    k = np.squeeze(mat['kspacelocations'])
    dcf = np.squeeze(mat['dcf'])
    kx = k.real 
    ky = k.imag 
    kerneltable = calc_kb_kernel(kwidth, overgrid_factor, klength)[0]
    nkernelpoints = kerneltable.size
    nsamples = k.size

    plt.plot(np.linspace(0, dcf.size, dcf.size), dcf)
    plt.show()

    grid.calcdcflut(kx, ky, dcf, kerneltable, gridsize, nsamples, nkernelpoints, kwidth)  

    plt.plot(np.linspace(0, dcf.size, dcf.size), dcf)
    plt.show()


