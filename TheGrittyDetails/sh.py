import numpy as np
import random
import math

class SHSample:
    def __init__(self):
        self.sph = None
        self.vec = None
        self.coeffs = None

def SH_setup_spherical_samples(sqrt_n_samples: int, n_bands: int):
    samples = []
    one_over_N = 1.0 / sqrt_n_samples

    counts = sqrt_n_samples*sqrt_n_samples
    counts_i = 0
    print('loops:', counts)
    for a in range(sqrt_n_samples):
        for b in range(sqrt_n_samples):
            # counter
            if counts_i % 100 == 0:
                print(' {}'.format(counts_i), end='')
                if counts_i > 0 and counts_i % 1000 == 0:
                    print()
            counts_i += 1

            sh = SHSample()

            x = (a + random.uniform(0, 1)) * one_over_N
            y = (b + random.uniform(0, 1)) * one_over_N
            theta = 2.0*math.acos(math.sqrt(1.0 - x))
            phi = 2.0 * np.pi * y
            r = 1.0
            sh.sph = np.float32([theta, phi, r])

            # spherical coords to unit vector
            sh.vec = np.float32([math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)])

            # precompute SH coeffecients
            max_sh_cnt = (n_bands-1)*((n_bands-1) + 1) + n_bands
            sh.coeffs = np.zeros(max_sh_cnt).astype(np.float32)
            for l in range(n_bands):
                for m in np.arange(-l, l+1, 1):
                    i = l*(l+1)+m
                    sh.coeffs[i] = SH(l, m, theta, phi)
            samples.append(sh)
    return samples

def P(l: int, m: int, x):
    # evaluate Associated Legendre Polynomial at x
    pmm = 1.0
    if m > 0:
        somx2 = math.sqrt((1.0 - x)*(1.0 + x))
        fact = 1.0
        for _ in range(1, m+1):
            pmm *= (-fact * somx2)
            fact += 2
    if l == m:
        return pmm

    pmmp1 = x * (2.0*m + 1.0) * pmm
    if l == (m+1):
        return pmmp1

    pll = 0.0
    for ll in range(m+2, l+1):
        pll = ((2.0*ll-1.0)*x*pmmp1-(ll+m-1.0)*pmm) / (ll-m)
        pmm = pmmp1
        pmmp1 = pll
    return pll

def K(l: int, m: int):
    # renormalization coefficient for SH function
    temp = ((2.0*l+1.0)*math.factorial(l-abs(m))) / (4.0*np.pi*math.factorial(l+abs(m)))
    return math.sqrt(temp)

def SH(l: int, m: int, theta, phi):
    """
    return a point sample of a SH basis function: y^m_l(theta, phi)
        l: band
        m: m < |l| int
        theta: [0, PI]
        phi: [0, 2PI]
    """
    if m == 0:
        return K(l, 0) * P(l, 0, math.cos(theta))
    elif m > 0:
        return math.sqrt(2.0) * K(l, m) * math.cos(m*phi) * P(l, m, math.cos(theta))
    else:
        return math.sqrt(2.0) * K(l, m) * math.sin(-m*phi) * P(l, -m, math.cos(theta))

def SH_project_polar_function(light_polar_fn, samples):
    """
    SH projection
    c_i = 4pi/N * sum_{j=1}^{N} light(x_j)*y_i(x_j)
    """
    n_coeffs = len(samples[0].coeffs)
    c = np.empty(n_coeffs).astype(np.float32)

    for sample in samples:
        theta = sample.sph[0]
        phi = sample.sph[1]
        L = light_polar_fn(theta, phi)
        for n in range(n_coeffs):
            c[n] += L * sample.coeffs[n]

    # divide the result by weight and number of samples
    weight = 4.0*np.pi
    factor = weight / len(samples)
    c = factor*c

    return c
