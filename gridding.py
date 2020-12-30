from scipy.io import loadmat
from kb_kernel import calc_kb_kernel
import matplotlib.pyplot as plt
import numpy as np
import grid

def read_klocations(filename):
    """ Reads a .txt containing k-space locations of data samples

        K-space must be normalized to be between [-0.5,0.5] 
        This is due to how the C-level implementation of the 
        gridding operation is done.
    """
    kx_locs = []
    ky_locs = []
    k = []
    with open(filename, 'r') as f:
        for line in f:
            kx_ky_tokens = line.split(',')
            kx_tokens = kx_ky_tokens[0].split('=')
            ky_tokens = kx_ky_tokens[1].split('=')

            k.append([float(kx_tokens[1].strip()), float(ky_tokens[1].strip())])

    k = np.asarray(k)

    return k

def normalize_klocations(k, f):
    """ Normalizes kx, ky values to be between [-f, f]
        
        f - normalize to this value, in most cases 1.
    """
    k[:,0] = np.multiply(k[:,0], f/np.nanmax(k[:,0]))
    k[:,1] = np.multiply(k[:,1], f/np.nanmax(k[:,1]))
    return k 

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

def combine_coil_data(im):
    """
    Combines the reconstructed image from all coil channels using
    the square root of the sum of the squares of the pixel values.
    """
    return np.sqrt(np.sum(np.square(im), axis=2))

if __name__ == "__main__":
    
    #mat = loadmat("data/spiralexampledata.mat")
    ball_phantom_mat = loadmat("./data/ball_phantom_data/ball_phantom.mat")

    k = read_klocations("./data/ball_phantom_data/rosetteKs.txt")
    k = normalize_klocations(k, 0.5)
    kx = k[:,0]
    ky = k[:,1]
    
    overgrid_factor = 8
    niterations = 1
    kwidth = 4.5 
    klength = 500
    gridsize = 256

    # squeeze mat arrays to remove extra dimensions
    #k = np.squeeze(mat['kspacelocations'])
    #dcf = np.squeeze(mat['dcf'])
    #s = np.squeeze(mat['spiraldata'])
    s = np.squeeze(ball_phantom_mat['data'])
    s = trim_kpoints(s)
    s = avg_kpoints(s)

    # Store reconstructed images for each coil
    ims = np.zeros((gridsize,gridsize,s.shape[1]), dtype=complex)

    kerneltable = calc_kb_kernel(kwidth, overgrid_factor, klength)[0]
    nkernelpoints = kerneltable.size
    nsamples = len(kx)

    dcf = np.zeros_like(kx)
    already_calculated = 0

    for i in range(s.shape[1]):
        s_temp = s[:,i,0]
        # Need to be in col major order (Fortran) to match with the way C code is written
        # This is not intuitive, since we would expect C major order.
        s_real = s_temp.real.flatten('F')
        s_imag = s_temp.imag.flatten('F')
        sg_real = np.zeros(gridsize**2)
        sg_imag = np.zeros(gridsize**2)

        if i > 0:
            already_calculated = 1
        #Calculate density compensation factors (DCF)
        dcf  = grid.calcdcflut(kx, ky, dcf, kerneltable, gridsize, nsamples, 
                                nkernelpoints, niterations, already_calculated, kwidth)  

        # Grid signal data
        s_real, s_imag = grid.gridlut(kx, ky, s_real, s_imag, dcf, sg_real, sg_imag, kerneltable, nsamples, gridsize, nkernelpoints, kwidth)

        sg = [complex(sg_real[i], sg_imag[i]) for i in range(sg_real.size)]
        sg = np.asarray(sg).reshape(gridsize,gridsize, order='F') 
        im = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(sg)))

        ims[:,:,i] = im
   
    im = combine_coil_data(ims)

    plt.imshow(abs(im), cmap='gray')
    plt.show()
