#!/usr/bin/python
from sys import stderr
from random import random, gauss, uniform
from math import pi, exp, sqrt, log
from subprocess import Popen, PIPE
from os import rename
from sys import argv

def warn(msg):
    stderr.write(msg + '\n')

# After Beaumont 2010, p. 388

# SMC ABC Parameters
N = 12500                                            # Number of samples
theta = list()                                      # Sets of parameter values that worked
omega = list()                                      # Weights associated with each parameter set
tau_sq = list()

prefix = ['./epi_sim_abc', '20', '10000']  # 
warn(' '.join(['#'] + prefix))

pars = ['R0', 'Ih', 'P0', 'h']

# Valid parameter ranges
R0_min = 1
R0_max = 5
Ih_min = 1.0/12.0 
Ih_max = 100
P0_min = 0 # Exponential, e.g. 2**P0_min
P0_max = 8 # 2**P0_max
h_min  = 0
h_max  = 1
h_mu   = 0.45
h_sigma= 0.2

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
    

def last_set_number():
    from glob import glob
    prior_files = glob('predictive_prior.*')
    prior_files = [ int(f.split('.')[1]) for f in prior_files if f.isdigit() ]
    if len(prior_files) == 0:
        return '0'
    else:
        # return max extension + 1, e.g. predictive_prior.3 --> 4
        return str(max(prior_files) + 1)


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
    for j in range(len(omega[-2])):
    #for j in range(N):
        denominator += omega[-2][j] \
                       * gauss_pdf(theta[-2]['R0'][j], sqrt(tau_sq[-1]['R0']), theta[-1]['R0'][-1]) \
                       * gauss_pdf(theta[-2]['Ih'][j], sqrt(tau_sq[-1]['Ih']), theta[-1]['Ih'][-1]) \
                       * gauss_pdf(theta[-2]['P0'][j], sqrt(tau_sq[-1]['P0']), theta[-1]['P0'][-1]) \
                       * gauss_pdf(theta[-2]['h'][j], sqrt(tau_sq[-1]['h']), theta[-1]['h'][-1]) 
    
    return numerator / denominator


# Example epi sim command:
# ./epi_sim_abc 100 10000 1.8 .3 25 .43 france.tab

# Read in last set of particles
theta.append( {'R0':[], 'Ih':[], 'P0':[], 'h':[] } )

for line in file("predictive_prior.1"):
    R0, Ih, h, P0 = [float(i) for i in line.split()]
    theta[0]['R0'].append( R0 )
    theta[0]['Ih'].append( Ih )
    theta[0]['P0'].append( P0 )
    theta[0]['h'].append( h )
 
omega.append( list() )
for w in file("predictive_prior_weights.1"):
    omega[0].append( float(w.strip()) ) 

last_num = last_set_number()
#rename("predictive_prior.out", "predictive_prior." + last_num)
#rename("weights.out", "weights." + last_num)

# Generate next set of particles
theta.append( {'R0':[], 'Ih':[], 'P0':[], 'h':[] } )
omega.append( list() )

tau_sq.append( dict() )
for par in ['R0', 'Ih', 'h', 'P0']:
    tau_sq[0][par] = 2 * var(theta[0][par]) # We will sample from kernels w/ 2x the variance of previous good params

output_tag = ''
if len(argv) > 1:
    output_tag = argv[1]
else:
    output_tag = str(int(last_num) + 1)

fo_weights = open('weights.' + output_tag, 'w')
fo_particles = open('particles.' + output_tag, 'w')
line = ','.join(['R0', 'Ih', 'h', 'P0', 'mean', 'median', 'sd', 'skew', 'ss', 'll', 'sl', 'ls', 'Re'])
fo_particles.write( line + '\n' )


for i in range(N):
    warn(output_tag + ' ' + str(i))
    idx = rand_nonuniform_int(omega[0])
    R0_mean = theta[0]['R0'][idx]
    Ih_mean = theta[0]['Ih'][idx]
    P0_mean = theta[0]['P0'][idx]
    h_mean  = theta[0]['h'][idx]

    R0 = trunc_gauss(R0_mean, sqrt(tau_sq[0]["R0"]), R0_min, R0_max)
    Ih = trunc_gauss(Ih_mean, sqrt(tau_sq[0]["Ih"]), Ih_min, Ih_max)
    P0 = round( 2**trunc_gauss(log(P0_mean, 2), log(sqrt(tau_sq[0]["P0"]), 2), P0_min, P0_max) )
    h  = trunc_gauss(h_mean, sqrt(tau_sq[0]["h"]), h_min, h_max)

    metrics_str  = Popen(prefix + [str(R0), str(Ih), str(P0), str(h), "france.tab"], stdout=PIPE).communicate()[0]
    metrics_iter = iter(metrics_str.split())
    metrics = dict(zip(metrics_iter, metrics_iter))
    
    line = ','.join([ metrics[k] for k in ['R0', 'Ih', 'h', 'P0', 'mean', 'median', 'sd', 'skew', 'ss', 'll', 'sl', 'ls', 'Re'] ])
    fo_particles.write( line + '\n' )

    theta[1]['R0'].append( R0 )
    theta[1]['Ih'].append( Ih )
    theta[1]['P0'].append( P0 )
    theta[1]['h'].append( h )

    particle_weight = weight(theta, omega, tau_sq)
    omega[1].append( particle_weight )
    fo_weights.write( str(particle_weight) + '\n' )

fo_weights.close()
fo_particles.close()
