INCLUDE = -I /home/tjhladish/work/EpiFire/src
LDFLAGS = /home/tjhladish/work/EpiFire/src/*.o

all: epi_sim_abc

epi_sim_abc: epi_sim_abc.cpp ../MultiSeason_Network/MultiSeason_Sim.h
	g++ -g -Wall $(INCLUDE) $(LDFLAGS) ../MultiSeason_Network/MultiSeason_Sim.h epi_sim_abc.cpp -o epi_sim_abc

