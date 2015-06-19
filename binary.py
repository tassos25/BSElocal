# module binary

from __future__ import print_function
from astropy import units as u
from astropy import constants as const
import numpy as np
import math
import timeit
import zams



def roche_lobe(m1, m2):
    if not isinstance(m1, u.Quantity):
        raise ValueError("m1 must be a Quantity with mass units")
    if not isinstance(m2, u.Quantity):
        raise ValueError("m2 must be a Quantity with mass units")

    q_mass1 = m1 / m2
    q_mass3 = q_mass1 ** (1.0 / 3.0)
    q_mass2 = q_mass1 ** (2.0 / 3.0)
    lobe = (0.49 * q_mass2) / (0.6 * q_mass2 + np.log(1.0 + q_mass3))
    return lobe  # Dimensionless. In units of orbital separation


def separation_to_period(separation, m1, m2):
    if not isinstance(m1, u.Quantity):
        raise ValueError("m1 must be a Quantity with mass units")
    if not isinstance(m2, u.Quantity):
        raise ValueError("m2 must be a Quantity with mass units")
    if not isinstance(separation, u.Quantity):
        raise ValueError("sep must be a Quantity with length units")

    mbin = m1 + m2
    period = np.sqrt(4.0 * math.pi ** 2.0 * separation ** 3.0 / (const.G * mbin))
    return period.to('day')


def period_to_separation(period, m1, m2):
    if not isinstance(m1, u.Quantity):
        raise ValueError("m1 must be a Quantity with mass units")
    if not isinstance(m2, u.Quantity):
        raise ValueError("m2 must be a Quantity with mass units")
    if not isinstance(period, u.Quantity):
        raise ValueError("period must be a Quantity with time units")

    mbin = m1 + m2
    separation = (period**2 * const.G * mbin / (4.0 * math.pi**2)) ** (1.0 / 3.0)
    return separation.to('Rsun')



def random_mass_ratio(nsample=1):
    return np.random.random(nsample)


def random_separation(amin, amax, nsample=1):
    """Returns a random number for the radius of separation
    of a binary. min cannot be zero
    """
    if not isinstance(amin, u.Quantity):
        raise ValueError("amin must be a Quantity with mass units")
    if not isinstance(amax, u.Quantity):
        raise ValueError("amax must be a Quantity with mass units")

    amin = amin.to('Rsun')
    amax = amax.to('Rsun')

    if amin <= 0.*u.Rsun:
        raise ValueError("min cannot be 0 or negative")
    if amax < amin:
        raise ValueError("min must be less than or equal to max")

    y =  np.random.random(nsample)*(np.log(amax.value) - np.log(amin.value))
    return np.exp(y)*amin





def random_eccentricity(nsample=1):
    # Generate a random number from the thermal eccentricity distribution f(e)=2e (CDF: F(e)=e^2)
    e = np.random.random(nsample)
    return np.sqrt(e)





def random_mass(pdf, pdf_comp, mmin, mmax, nsample=1):
    """Generates a random value for the mass of a star
    where min and max define the range for the value of
    the mass and a is the exponent associated with the
    probability distribution for masses between mmin and mmax
    solar masses. mmin can be no smaller than 0.08 and mmax
    can be no more than 150. a must be within the range
    2.35-3.2."""

    if not isinstance(mmin, u.Quantity):
        raise ValueError("mmin must be a Quantity with mass units")
    if not isinstance(mmax, u.Quantity):
        raise ValueError("mmax must be a Quantity with mass units")

    mmin = mmin.to('Msun')
    mmax = mmax.to('Msun')




    # Counters
    naccept = 0
    ntrial = 0

    # Keeps generating numbers until we achieve the desired n
    ran=[] # output list of random numbers
    
    c = pdf(mmin.value)/pdf_comp(mmin.value)
    
    while naccept < nsample:
        uu = np.random.random()
        v = random_power(min=mmin.value, max=mmax.value, slope=-1.3)
        ntrial += 1
        while c*uu > f(v)/g(v):
            # reject v and pick a new one
            uu = np.random.random()
            v = random_power(min=mmin.value, max=mmax.value, slope=-1.3)
            ntrial += 1
        ran.append(v)
        naccept += 1

    ran = np.asarray(ran)

    return ran * u.Msun




def f(x, a1 = 1.3, a2 = 2.3, a3 = 2.35, m_min=0.08, m_1 = 0.5, m_2 = 1.0, m_max=150.0):
    fo_inv = 1./(1.-a1) * (m_1**(1.-a1)-m_min**(1.-a1))/(m_min**(-a1)) \
        + 1./(1.-a2) * (m_1/m_min)**(-a1) * (m_2**(1.-a2)-m_1**(1.-a2))/(m_1**(-a2)) \
        + 1./(1.-a3) * (m_1/m_min)**(-a1) * (m_2/m_1)**(-a2) * (m_max**(1.-a3)-m_2**(1.-a3))/(m_2**(-a3))

    fo = 1. / fo_inv

    if type(x) is np.ndarray:
        value = np.zeros(len(x))
        idx1 = np.where((x >= m_min) * (x < m_1))
        if len(idx1[0]) > 0:
            value[idx1] = fo * np.power(x[idx1]/m_min,-a1)
        idx2 = np.where((x >= m_1)  * (x < m_2))
        if len(idx2[0]) > 0:
            value[idx2] = fo * np.power(m_1/m_min,-a1) * np.power(x[idx2]/m_1,-a2)
        idx3 = np.where((x >= m_2) * (x < m_max))
        if len(idx3[0]) > 0:
            value[idx3] = fo * np.power(m_1/m_min,-a1) * np.power(m_2/m_1,-a2) * np.power(x[idx3]/m_2,-a3)
        idx4 = np.where(x > m_max)
        if len(idx4[0]) > 0:
            value[idx4] = 0.0
        idx5 = np.where(x < m_min)
        if len(idx5[0]) > 0:
            value[idx5] = 0.0

    elif type(x) is float or type(x) is np.float64:
        value=0.0
        if (x >= m_min) and (x < m_1):
            value = fo * np.power(x/m_min,-a1)
        elif (x >= m_1)  and (x < m_2):
            value = fo * np.power(m_1/m_min,-a1) * np.power(x/m_1,-a2)
        elif (x >= m_2) and (x < m_max):
            value = fo * np.power(m_1/m_min,-a1) * np.power(m_2/m_1,-a2) * np.power(x/m_2,-a3)
        else:
            value=0.0

    return value



def g(x, a1 = 1.3, m_min=0.08, m_max=150.0):
    fo_inv = 1./(1.-a1) * (m_max**(1.-a1)-m_min**(1.-a1))/(m_min**(-a1)) 
 
    fo = 1. / fo_inv

    if type(x) is np.ndarray:
        value = np.zeros(len(x))
        idx1 = np.where((x >= m_min) * (x < m_max))
        if len(idx1[0]) > 0:
            value[idx1] = fo * np.power(x[idx1]/m_min,-a1)
    elif type(x) is float or type(x) is np.float64:
        value=0.0
        if (x >= m_min) and (x < m_max):
            value = fo * np.power(x/m_min,-a1)

    return value


def random_power(min, max, slope):
    """Returns a random number with a power law distribution of power n,      
    from values between min and max."""

    range = np.power(max, 1+slope) - np.power(min, 1+slope)
    return np.power(np.random.random()*range + np.power(min, 1+slope), 1.0/(1+slope))




def generate_bse_binaries(fileout, nbin):
    time_to_evolve = 100. * u.Myr  #Myr
    z=0.006
    zams.prepare_coefficients(z)
    sep_max = 1.e5 * u.Rsun  #Rsun
    fout = open(fileout,'w')
    print(nbin, file=fout)
    
    for i in range(nbin):
        m1 = random_mass(f, g, 15.0*u.Msun, 100.*u.Msun)
        q_ratio = np.random.random()
        m2 = m1 * q_ratio
        ecc = random_eccentricity()
        rl1 = roche_lobe(m1, m2)
        rzams1 = zams.rzams(m1)
        sep_min = rzams1/(rl1*(1.-ecc))
        # Sometime eccentricity gets too close to 1 (eg. 0.999997) and no binary can be formed
        # Instead of limiting ecc to a max value, we just make another draw
        while sep_min > sep_max:
            ecc = random_eccentricity()
            sep_min = rzams1/(rl1*(1.-ecc))
            
        sep = random_separation(sep_min, sep_max)
        tb = separation_to_period(sep, m1, m2)
        print(float(m1.value), float(m2.value), float(tb.value), float(ecc), z, time_to_evolve.value, file=fout)
        
    fout.close()




if __name__ == "__main__":
	generate_bse_binaries('binariesBH1e7.in', 10000000)
