EPIFIRE = "/home/tjhladish/work/EpiFire/src/"

msms_sim: msms_sim.cpp MultiSeason_MultiStrain_Sim.h
	g++ -O2 -I $(EPIFIRE)  $(EPIFIRE)*.o MultiSeason_MultiStrain_Sim.h msms_sim.cpp -o msms_sim

msmsc_sim: msmsc_sim.cpp MSMSC.h
	g++ -O2 -I $(EPIFIRE)  $(EPIFIRE)*.o MSMSC.h msmsc_sim.cpp -o msmsc_sim

