CC=g++
CFLAGS = -O2
INCLUDE = -I $(HOME)/EpiFire/src
LDFLAGS = $(HOME)/EpiFire/src/*.o

all: epi_sim_abc

epi_sim_abc: epi_sim_abc.cpp ../MultiSeason_Network/MultiSeason_Sim.h
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) epi_sim_abc.cpp -o epi_sim_abc

