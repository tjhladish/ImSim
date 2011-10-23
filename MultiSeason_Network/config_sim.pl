#CONFIG FILE 

use vars qw/ $sim_exec $j_max $reuse_net $s_max $n $R0 $dist $p1 $p2 $Ih $model $m $P0_ct $P0_model $RunID $EdgeFile/;

print "# Reading epidemic simulation config file (config_sim.pl) ...\n";

#our $path      = '/home/tjhladish/work/EpiFire';
our $sim_exec  = '/home/tjhladish/work/ImSim/MultiSeason_Network/epi_sim';

our $j_max     = 5000;     # number of different networks to iterate through ***Or reps when reuse_net == true***
our $reuse_net = 1;        # Do you reuse the network or generate a new one for each j?
our $s_max     = 50;       # number of seasons to simulate
our $n         = 25622;    # size of network to use
our $R0        = 2.4;      # R value for naive network

our $dist      = 3;        # Dist: 0 = poisson, 1 = exponential, 2 = powerlaw, 3 = urban, 4 = constant, 5 = use edgelist file
our $p1        = -1;       # Distribution parameter 1 (lamda for poisson, exponential; something else for powerlaw)
our $p2        = -1;       # Dist. parameter 2 (used for powerlaw)
our $Ih        = 5;        # Immunity halflife, in years
our $model     = 0;        # Mixture Model: 0 = none, 1 = half, 2 = full, 3 = natural, 4 = half_shuffle, 5 = full_shuffle
our $m         = 0.1;      # Fraction of individuals that migrate in natural mixture model

our $P0_ct     = 10;       # Number of patient zeros/infections to introduce (usually 1)
our $P0_model  = 0;        # What patient-zero model? 0 = static number of randomly drawn nodes, 1 = dynamic number (sensitive 
                           #   to node immunity) of randomly drawn nodes, 2 = static number of clustered nodes
our $RunID     = 7;        # Arbitrary ID for bookkeeping purposes--you might want to keep it short, though
our $EdgeFile  = "";#10000_household_undirected.edges";

print "# sim_exec  = $sim_exec\n";
print "# j_max     = $j_max\n";
print "# reuse_net = $reuse_net\n";
print "# s_max     = $s_max\n";
print "# n         = $n\n";
print "# R0        = $R0\n";
print "# dist      = $dist\n";
print "# p1        = $p1\n";
print "# p2        = $p2\n";
print "# Ih        = $Ih\n";
print "# model     = $model\n";
print "# m         = $m\n";
print "# P0_ct     = $P0_ct\n";
print "# P0_model  = $P0_model\n";
print "# RunID     = $RunID\n";
print "# EdgeFile  = $EdgeFile\n";
print "# ... done.\n";
