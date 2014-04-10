#!/usr/bin/python
from numpy import mean, median
from random import random

def mean_max(d):
    vals = []
    for loc in d.keys():
        vals.append(max(d[loc].values()))
    return mean(vals)

def mean_range(d):
    vals = []
    for loc in d.keys():
        vals.append(max(d[loc].values()) - min(d[loc].values()))
    return mean(vals)


def import_observed_data(filename):
    obs = dict()
    obs_flattened = []
    for line in file(filename):
        loc, year, val = line.strip().split(',')
        obs_flattened.append(float(val))
        if loc not in obs:
            obs[loc] = dict()
        obs[loc][int(year)] = float(val)
    return obs, obs_flattened

def skew(d):
    return mean(d) - median(d)

def stdev(d):
    x = mean(d)
    var = sum( [(x - y)**2 for y in d] ) / (len(d) - 1)
    return var**0.5

def autocorrelation_matrix(data, randomize_start=False):
    ac_matrix = {'ss':0, 'll':0, 'sl':0, 'ls':0}
    N = 0.0

    for loc in data.keys():
        if len(data[loc]) <= 1: continue

        avg = mean(data[loc].values())
        sorted_years = sorted(data[loc].keys()[:-1]) 

        if randomize_start == True:
            if random() > 0.5:
                sorted_years = sorted_years[:-1]
            if random() > 0.5:
                sorted_years = sorted_years[1:]

        for year in sorted_years:
            this_year = data[loc][year]
            next_year = data[loc][year+1]
            if this_year < avg and next_year < avg:
                ac_matrix['ss'] += 1.0
            elif this_year >= avg and next_year >= avg:
                ac_matrix['ll'] += 1.0 
            elif this_year < avg and next_year >= avg:
                ac_matrix['sl'] += 1.0
            else:
                ac_matrix['ls'] += 1.0
            N += 1.0
    for k in ac_matrix.keys():
        ac_matrix[k] = ac_matrix[k]/N
    return ac_matrix

obs, obs_flattened = import_observed_data('paris.csv')
ac = autocorrelation_matrix(obs, False)
print "mean, median, max, range, sd, skew, ac_ss, ac_ll, ac_sl, ac_ls"
print mean(obs_flattened), median(obs_flattened), mean_max(obs), mean_range(obs), stdev(obs_flattened), skew(obs_flattened), ac['ss'], ac['ll'], ac['sl'], ac['ls'],
#for i in range(1000):
#    ac = autocorrelation_matrix(obs, True)
#    print ac['ss'], ac['ll'], ac['sl'], ac['ls']
