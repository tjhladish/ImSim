#!/usr/bin/python
from numpy import mean, median

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

def autocorrelation_matrix(data):

    N = 0.0;

    # Determine total number of comparisons available
    for k in data:
        if len(data[k]) > 1:
            N += len(data[k]) - 1    # n-1 per time series

    Inv_N = 1.0 / float(N)          # percent contributed by each comparison

    # Autocorrelation scores
    # indexing: 0 == small, small; 1 == large, large; 2 == small, large; 3 == large, small
    ac_matrix = {'ss':0, 'll':0, 'sl':0, 'ls':0}

    for loc in data.keys():
        if len(data[loc]) <= 1: continue

        avg = mean(data[loc].values())
        sorted_years = sorted(data[loc].keys()[:-1]) 
        for year in sorted_years:
            this_year = data[loc][year]
            next_year = data[loc][year+1]
            if this_year < avg and next_year < avg:
                ac_matrix['ss'] += Inv_N
            elif this_year >= avg and next_year >= avg:
                ac_matrix['ll'] += Inv_N
            elif this_year < avg and next_year >= avg:
                ac_matrix['sl'] += Inv_N
            else:
                ac_matrix['ls'] += Inv_N
    return ac_matrix;

obs, obs_flattened = import_observed_data('france.tab')
ac = autocorrelation_matrix(obs)
print "mean, median, sd, skew, ac_ss, ac_ll, ac_sl, ac_ls"
print mean(obs_flattened), median(obs_flattened), stdev(obs_flattened), skew(obs_flattened), ac['ss'], ac['ll'], ac['sl'], ac['ls'],
