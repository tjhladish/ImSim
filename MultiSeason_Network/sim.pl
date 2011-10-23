#!/usr/bin/perl  -w
use strict;

BEGIN { require "config_sim.pl"; }

my $cmd = "$sim_exec $j_max $reuse_net $s_max $n $R0 $dist $p1 $p2 $Ih $model $m $P0_ct $P0_model $RunID $EdgeFile";
warn "# $cmd\n";
system($cmd);
