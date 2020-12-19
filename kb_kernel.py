import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv


def kb(u, w, beta):
    """
    Computes the Kaiser-Bessel function used for gridding, namely:

    Based on Jackson et al.

    y = f(u,w,beta) = I_0 [beta * sqrt(1-(2u/w)^2)]/w

    where I_0 is the zero-order modified Bessel function of the first kind.

    Params:
        u - vector of k-space locations for calculation
        w - width parameter
        beta - beta parameter
    """

    y = np.zeros(u.shape)
    # iv is the modified Bessel function
    x = lambda u:  iv(0, beta * np.sqrt(1 - (2*u/w)**2)) / w
    # Only calc. values within kernel width; otherwise, set to zero. 
    y = [x(ui) if ui < w/2 else 0 for ui in u]

    return y
    

def calc_kb_kernel(kwidth, overgrid_factor, klength):
    """
    Calculates the appropriate Kaiser-Bessel kernel for gridding, using the approach
    of Jackson et al. 

    Params:
        kwidth - kernel width in grid samples
        overgrid_factor - the factor by which the grid oversamples the original acquisition grid.
        klength - kernel look-up-table length. 

    Returns:
        kern - kernel values for klength values of u, uniformly spaced from 1 to kwidth/2.
        u - u values.
    """
    
    a = overgrid_factor
    w = kwidth
    beta = np.pi * np.sqrt(w**2 / a**2 * (a-0.5)**2 - 0.8)

    u = np.linspace(0, klength, klength) * (kwidth/(2*klength))

    kern = kb(u, kwidth, beta)
    kern = kern / max(kern) # normalize

    return kern, u

if __name__ == "__main__":

    kernel, u = calc_kb_kernel(1.5, 2, 1000)
    plt.scatter(u, kernel)
    plt.show()

