#!/usr/bin/python
from sys import stderr, argv, exit
from random import random, gauss, uniform, randint
from math import pi, exp, sqrt, log
from subprocess import Popen, PIPE

def warn(msg):
    stderr.write(msg + '\n')

# After Beaumont 2010, p. 388

# Valid parameter ranges
R0_min = 1
R0_max = 8
Ih_min = 1.0/12.0 
Ih_max = 100
P0_min = 1
P0_max = 256
h_min  = 0
h_max  = 1


def read_config():
    config = dict()
    for line in file('abc.config'):
        k, v = line.strip().split('\t')
        try:
            v = int(v)
        except ValueError:
            try:
                v = float(v)
            except ValueError:
                pass
        config[k] = v
    return config

def valid_pars(R0, Ih, P0, h):
    if R0 < R0_min or R0 > R0_max:
        return false
    if Ih < Ih_min  or Ih > Ih_max:
        return false
    if P0 < P0_min or P0 > PO_max:
        return false 
    if h < h_min   or h > h_max:
        return false
    return true
    

def sample_priors():
    # R0, Ih, P0, h
    R0 = uniform(R0_min, R0_max)
    Ih = uniform(Ih_min, Ih_max) # infinite immunity <-> losing 99% of immunity every year
    h  = uniform(h_min, h_max)
    P0 = randint(P0_min, P0_max)
    #P0 = round(2**uniform(P0_min,P0_max)) # [1,256] loguniform
    #h  = trunc_gauss(h_mu, h_sigma, h_min, h_max)
    return [R0, Ih, h, P0]


def trunc_gauss(mu, sigma, a, b):
    while 1:
        dev = gauss(mu, sigma) # mean and sd
        if dev > a and dev < b:
            return dev

def var(s):
    v = 0
    x = sum(s)/float(len(s))
    for i in range(len(s)):
        v += (s[i] - x)**2
    return v

def rand_nonuniform_int(w):
    s = sum(w)
    x = random()
    r = w[0]/float(s)
    i = 0
    while r < x:
        i += 1
        r += w[i]/float(s)
    return i

def uniform_pdf(a,b):
    return 1.0/abs(b-a)

def gauss_pdf(mu, sigma, x):
    p = 1.0/sqrt(2*pi*sigma**2) * exp(-(x-mu)**2/ (2*sigma**2) )
    return p

def weight(theta, omega, tau_sq):
    numerator = uniform_pdf(R0_min, R0_max) \
                * uniform_pdf(Ih_min, Ih_max) \
                * uniform_pdf(P0_min, P0_max) \
                * gauss_pdf(h_mu, h_sigma, theta[-1]['h'][-1])

    denominator = 0
    for j in range(N):
        denominator += omega[-2][j] \
                       * gauss_pdf(theta[-2]['R0'][j], sqrt(tau_sq[-1]['R0']), theta[-1]['R0'][-1]) \
                       * gauss_pdf(theta[-2]['Ih'][j], sqrt(tau_sq[-1]['Ih']), theta[-1]['Ih'][-1]) \
                       * gauss_pdf(theta[-2]['P0'][j], sqrt(tau_sq[-1]['P0']), theta[-1]['P0'][-1]) \
                       * gauss_pdf(theta[-2]['h'][j], sqrt(tau_sq[-1]['h']), theta[-1]['h'][-1]) 
    
    return numerator / denominator

conf = read_config()
prefix = ['./epi_sim_abc', str(conf['burnin']), str(conf['network_size'])]
particle_keys = ['R0', 'Ih', 'h', 'P0', 'mean', 'median', 'max', 'range', 'sd', 'skew', 'ss', 'll', 'sl', 'ls', 'Re']

fo_parameters = open('parameters.0', 'w')

warn(' '.join(['#'] + prefix))
for i in range(conf['sample_size']):
    #warn(output_tag + ' ' + str(i))
    par_str = [str(i) for i in sample_priors()]
    line = ','.join(par_str)
    fo_parameters.write( line + '\n' )

fo_parameters.close()
