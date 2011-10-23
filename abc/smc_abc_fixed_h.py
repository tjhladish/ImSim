#!/usr/bin/python
import random 
import math
import subprocess
import pp

# After Beaumont 2010, p. 388

T = 1
N = 10000

epsilon_max = 10
epsilon_min = 1

epsilon = [epsilon_min + i*(epsilon_max-epsilon_min)/T for i in range(T,-1,-1)]
theta = list()
omega = list()
tau_sq = list()

prefix = ['./epi_sim_abc', '1', '1', '50', '10000']
pars = ['R0', 'D', 'P0']

# Valid parameter ranges
lim = { 'R0':
            { 'min': 1,
              'max': 5},
        'D':
            { 'min': 0,
              'max': 10}, # 6.64385619
        'P0':
            { 'min': 2**0, # Exponential, e.g. 2**P0_min
              'max': 2**6} # 2**P0_max
       }

#h = 0.45

# tuple of all parallel python servers to connect with
ppservers = ()

# Creates jobserver with automatically detected number of workers
job_server = pp.Server(ppservers=ppservers,secret="123")
job_server.set_ncpus(6)


def sample_priors():
    # R0, D, P0
    R0 = random.uniform(lim['R0']['min'], lim['R0']['max'])
    D  = random.uniform(lim['D']['min'], lim['D']['max']) # infinite immunity <-> losing 99% of immunity every year
    P0 = round(2**random.uniform( math.log(lim['P0']['min'],2), math.log(lim['P0']['max'],2) )) # [1,256] loguniform

    return R0, D, P0


def trunc_gauss(mu, sigma, par_range):
    while 1:
        dev = random.gauss(mu, sigma) # mean and sd
        if dev > par_range['min'] and dev < par_range['max']:
            return dev

def mean(s):
    return sum(s)/float(len(s))

def var(s):
    v = 0
    x = mean(s)
    for i in range(len(s)):
        v += (s[i] - x)**2
    return v

def rand_nonuniform_int(w):
    s = sum(w)
    x = random.random()
    r = w[0]/float(s)
    i = 0
    while r < x:
        i += 1
        r += w[i]/float(s)
    return i

def uniform_pdf(par_range):
    a,b = par_range['min'], par_range['max']
    return 1.0/abs(b-a)

def gauss_pdf(mu, sigma, x):
    p = 1.0/math.sqrt(2*math.pi*sigma**2) * math.exp(-(x-mu)**2/ (2*sigma**2) )
    return p

def weight(theta, omega, tau_sq):
    numerator = uniform_pdf(lim['R0']) \
                * uniform_pdf(lim['D']) \
                * uniform_pdf(lim['P0'])

    denominator = 0
    for j in range(N):
        denominator += omega[-2][j] \
                       * gauss_pdf(theta[-2]['R0'][j], math.sqrt(tau_sq[-1]['R0']), theta[-1]['R0'][-1]) \
                       * gauss_pdf(theta[-2]['D'][j], math.sqrt(tau_sq[-1]['D']), theta[-1]['D'][-1]) \
                       * gauss_pdf(theta[-2]['P0'][j], math.sqrt(tau_sq[-1]['P0']), theta[-1]['P0'][-1])
    
    return numerator / denominator
#http://www.parallelpython.com/content/view/17/31/

def run_sim(prefix, R0, D, P0, obs_file):
    h = 0.45
    d = subprocess.Popen(prefix + [str(R0), str(D), str(P0), str(h), "france.tab"], stdout=subprocess.PIPE).communicate()[0]
    d = float(d.strip())
    return d

def find_good_theta(t, epsilon, omega, theta, tau_sq, lim, prefix):
    KS = ()
    while KS >= epsilon[t]:
        idx = rand_nonuniform_int(omega[t-1])
        R0_mean = theta[t-1]['R0'][idx]
        D_mean = theta[t-1]['D'][idx]
        P0_mean = theta[t-1]['P0'][idx]

        R0 = trunc_gauss(R0_mean, math.sqrt(tau_sq[t-1]["R0"]), lim['R0'])
        D = trunc_gauss(D_mean, math.sqrt(tau_sq[t-1]["D"]), lim['D'])
        P0 = round(trunc_gauss(P0_mean, math.sqrt(tau_sq[t-1]["P0"]), lim['P0']))

        exec_str = ' '.join(prefix + [str(R0), str(D), str(P0), 'h', "france.tab"])
        #print exec_str
        KS = run_sim(prefix, R0, D, P0, "france.tab")
        #KS_all.append(KS)
#        print "*IN* KS, t, i, R0, D, P0:", '\t'.join([str(j) for j in [KS, t, i, R0, D, P0]])
    return KS, R0, D, P0

#  Example epi sim command:
# ./epi_sim_abc 1 1 100 10000 1.8 .3 25 .43 france.tab

# Step 1
theta.append( dict() )
for par in ['R0', 'D', 'P0']:
    theta[0][par] = []
omega.append( list() )

KS_vals = []
KS_all = []
print "\nKS threshold:", epsilon[0]
for i in range(N):
    KS = ()
    R0, D, P0 = 0,0,0
    while KS >= epsilon[0]:
        R0, D, P0 = sample_priors()
        exec_str = ' '.join(prefix + [str(R0), str(D), str(P0), 'h', "france.tab"])
        #print exec_str
        KS = run_sim(prefix, R0, D, P0, "france.tab")
        KS_all.append(KS)
    KS_vals.append(KS)
    print "KS, t, i, R0, D, P0:", '\t'.join([str(j) for j in [KS, 0, i, R0, D, P0]])
    
    theta[0]['R0'].append( R0 )
    theta[0]['D'].append( D )
    theta[0]['P0'].append( P0 )
    
    omega[0].append( 1.0/N )

#print "step 1, mean KS (all):", mean(KS_all), ", mean KS (best):", mean(KS_vals)
tau_sq.append(dict()) # initialize dictionary
for par in ['R0', 'D', 'P0']:
    tau_sq[0][par] = 2 * var(theta[0][par]) # We will sample from kernels w/ 2x the variance of previous good params


# Steps 2 through T
for t in range(1,T):
    print "\nKS threshold:", epsilon[t]
    KS_vals = []
    KS_all = []
    theta.append( dict() )
    for par in ['R0', 'D', 'P0']:
        theta[t][par] = []
    omega.append( list() )
    jobs = []
    for i in range(N):
        jobs.append( job_server.submit(find_good_theta,(t, epsilon, omega, theta, tau_sq, lim, prefix),\
                            (rand_nonuniform_int, trunc_gauss, run_sim),\
                            ("math", "random", "subprocess")) )

    i=0
    for job in jobs:
        KS, R0, D, P0 = job()
        #KS, R0, D, P0 = find_good_theta(t, epsilon, omega, theta, tau_sq, lim, prefix);
        KS_vals.append(KS)
        print "KS, t, i, R0, D, P0:", '\t'.join([str(j) for j in [KS, t, i, R0, D, P0]])
        i+=1

        theta[t]['R0'].append( R0 )
        theta[t]['D'].append( D )
        theta[t]['P0'].append( P0 )

        omega[t].append( weight(theta, omega, tau_sq) )
    
    #print "step", str(t) + ", mean KS (all):", mean(KS_all), ", mean KS (best):", mean(KS_vals)
    print theta[-1]
    tau_sq.append(dict()) # initialize dictionary
    for par in ['R0', 'D', 'P0']:
        tau_sq[-1][par] = 2 * var(theta[-1][par]) # We will sample from kernels w/ 2x the variance of previous good params

#########################################
for i in range(len(theta)):
    print str(i) + ':\t',
    for par in ['R0', 'D', 'P0']:
        print theta[i][par][-1],
    print '\n'


fo = open('best_par', 'w')
for par in ['R0', 'D', 'P0']:
    vals = '\t'.join([str(i) for i in theta[-1][par]])
    fo.write(par + '\n' + vals + '\n')









