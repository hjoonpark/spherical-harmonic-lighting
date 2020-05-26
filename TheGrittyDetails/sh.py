import numpy as np
import random
import math
import vtk

class SHSample:
    def __init__(self):
        self.sph = None
        self.vec = None
        self.coeffs = None

def load_obj(obj_path):
    with open(obj_path, 'r') as f:
        lines = f.readlines()
        vs = []
        vts = []
        vns = []
        faces = []
        f_idx = 0
        for line in lines:
            l = line.split(' ')
            if l[0] == 'v':
                vs.append([float(l[1]), float(l[2]), float(l[3])])
            elif l[0] == 'vt':
                vts.append([float(l[1]), float(l[2])])
            elif l[0] == 'vn':
                vns.append([float(l[1]), float(l[2]), float(l[3])])
            elif l[0] == 'f':
                fis = []
                for i in range(1, len(l)):
                    v = l[i].split('/')
                    fi = [int(v[0]), int(v[1]), int(v[2])]
                    fis.append(fi)
                faces.append(fis)
                f_idx += 1
        f.close()
    return {'vs': np.float32(vs), 'vns': np.float32(vns), 'vts': np.float32(vts), 'faces': faces}

def get_unique_vertex_normal_pairs(obj_data):
    """
    f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3

    out: vertex_data = list of dict: [{'vertex': 3d, 'normal': 3d}]
    """
    faces = obj_data['faces']
    vs = obj_data['vs']
    vns = obj_data['vns']
    vns_avged = np.zeros(vs.shape).astype(np.float32)
    vns_count = np.zeros(len(vs)).astype(np.uint8)

    for f in faces:
        for fi in f:
            vertex_idx = fi[0] - 1
            normal_idx = fi[2] - 1
            vns_avged[vertex_idx] += vns[normal_idx]
            vns_count[vertex_idx] += 1

    for i in range(vns_avged.shape[0]):
        vns_avged[i] /= vns_count[i]

    vertex_data = []
    for i in range(vs.shape[0]):
        v = {'vertex': vs[i], 'normal': vns_avged[i]}
        vertex_data.append(v)
    return vertex_data

def compute_transfer_matrix(samples, n_coeff: int, vertex_data):
    """
    vertex_data = [{'vertex': 3d, 'normal': 3d}]
    """
    albedo = np.float32([1, 1, 1]) # assume white texture for the obj's surface

    for v_idx, v_data in enumerate(vertex_data):
        normal = v_data['normal']
        coeff = np.zeros((n_coeff, 3)).astype(np.float32) # 3 for rgb

        for sample in samples: # sum over all samples of light source
            H = sample.vec.dot(normal)
            if H > 0.0:
                # light ray inside hemisphere centered at the vertex point. SH project over all bands into the sum vector
                for j in range(n_coeff):
                    c_j = H * sample.coeffs[j]
                    coeff[j, 0] += albedo[0] * c_j
                    coeff[j, 1] += albedo[1] * c_j
                    coeff[j, 2] += albedo[2] * c_j
        # device the result by probability / number of samples
        factor = (4.0*np.pi) / len(samples)
        for j in range(n_coeff):
            coeff[j, :] *= factor
        vertex_data[v_idx]['transfer_vector'] = coeff

    return vertex_data




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
        return math.sqrt(2.0) * K(l, -m) * math.sin(-m*phi) * P(l, -m, math.cos(theta))

def SH_project_polar_function(polar_fn, samples):
    """
    SH projection
    c_i = 4pi/N * sum_{j=1}^{N} light(x_j)*y_i(x_j)
      i: SH coefficient index
      N: number of samples
    """
    n_coeffs = len(samples[0].coeffs)
    c = np.empty(n_coeffs).astype(np.float32)

    for sample in samples:
        theta = sample.sph[0]
        phi = sample.sph[1]
        L = polar_fn(theta, phi)
        for n in range(n_coeffs):
            c[n] += L * sample.coeffs[n]

    # divide the result by weight and number of samples
    weight = 4.0*np.pi
    factor = weight / len(samples)
    c = factor*c

    return c