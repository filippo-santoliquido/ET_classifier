import numpy as np

def uniform_pdf(x, a, b):
    return np.where(np.logical_and(x >= a, x <= b), 1 / (b - a), 0)

def uniform_in_cosine_pdf(x, a = -np.pi/2, b = np.pi/2):
    return np.where(np.logical_and(x >= a, x <= b), 0.5 * np.cos(x), 0)

def uniform_in_sine_pdf(x, a = 0, b = np.pi):
    return np.where(np.logical_and(x >= a, x <= b), 0.5 * np.sin(x), 0)

import astropy
from astropy.cosmology import Planck18 as cosmo
from astropy.cosmology import z_at_value
from astropy import units as u
from astropy import constants as const

"""
This is Ulyana's orginal version
def norm_factor_distance_prior(a, b):
    z_low = z_at_value(cosmo.luminosity_distance, a * u.Mpc).value
    z_high = z_at_value(cosmo.luminosity_distance, b * u.Mpc).value
    zz = np.linspace(z_low, z_high, 1000)
    xx = cosmo.luminosity_distance(zz).value
    yy = cosmo.differential_comoving_volume(zz).value * (xx / (1 + zz) + const.c.value * (1 + zz) / (cosmo.H(0).value * cosmo.efunc(zz)))**(-1)/ (1 + zz)
    return np.trapz(yy, zz)
"""

"""
This is Ulyana's orginal version
def uniform_in_differential_comoving_volume_pdf(x, a, b, z, norm_factor):
    #norm_factor = norm_factor_distance_prior(a, b)
    return np.where(np.logical_and(z >= a, z <= b), 
                cosmo.differential_comoving_volume(z).value * (x / (1 + z) + (const.c.value / 1000) * (1 + z) / (cosmo.H(0).value * cosmo.efunc(z)))**(-1) / (1 + z)/norm_factor, 0)
"""


##### This is my version - exaclty the same as above, just written differently and normalised - #####
def uniform_in_differential_comoving_volume_pdf(x, a, b, z, norm_factor):
    
    dVc_dz = 4 * np.pi * cosmo.differential_comoving_volume(z).value
    
    dl = x #cosmo.luminosity_distance(z).value
    
    dz_ddl =  ( (const.c.value / 1000) * ( (z+1)/(cosmo.H(0).value * cosmo.efunc(z)) + dl/(1+z)/(const.c.value / 1000)) )**(-1)
    
    dVc_ddl = dVc_dz * dz_ddl * (1/(1+z)) / norm_factor
    
    return np.where(np.logical_and(z >= a, z <= b), dVc_ddl , 0)

##### This is my version - probably correct - #####
def norm_factor_distance_prior(z_low, z_high):
    
    zz = np.geomspace(z_low, z_high, 1000)

    dVc_dz = 4 * np.pi * cosmo.differential_comoving_volume(zz).value
    
    dl = cosmo.luminosity_distance(zz).value
    
    dz_ddl =  ( (const.c.value / 1000) * ( (zz+1)/(cosmo.H(0).value * cosmo.efunc(zz)) + dl/(1+zz)/(const.c.value / 1000)) )**(-1)
    
    dVc_ddl = dVc_dz * dz_ddl * (1/(1+zz)) 
    
    return np.trapz(dVc_ddl, zz)


# Sample from priors

def uniform_sampler(lower_lim, upper_lim, N_points):
    return np.random.uniform(lower_lim, upper_lim, N_points)

def from_Mchirp_q_to_masses(Mchirp, q):
    m1 = Mchirp * (1 + q)**(1/5) * q**(3/5)
    m2 = Mchirp * (1 + q)**(1/5) * q**(-2/5)
    return m1, m2

def from_m1_m2_Mchirp_q(mass_1, mass_2):
    Mchirp = (mass_1 * mass_2)**(3/5) / (mass_1 + mass_2)**(1/5)
    q = mass_2 / mass_1
    return Mchirp, q

def chirp_mass(mass_1, mass_2):
    return (mass_1 * mass_2)**(3/5) / (mass_1 + mass_2)**(1/5)

def mass_sampler_uniform_in_q_and_chirp(mmin, mmax, N_points):
    chirp_low = chirp_mass(mmin, mmin)
    chirp_high = chirp_mass(mmax, mmax)
    Mchirp = uniform_sampler(chirp_low, chirp_high, N_points)
    q = uniform_sampler(0, 1, N_points)
    m1, m2 = from_Mchirp_q_to_masses(Mchirp, q)
    return m1, m2

def uniform_in_cosine_sampler(low, high, N_points):
    if low < -1 or high > 1:
        print('Error: bounds must be between -1 and 1')
    return np.arccos(np.random.uniform(low, high, N_points))

