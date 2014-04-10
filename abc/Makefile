CC=g++
CFLAGS = -O2
INCLUDE = -I $(HOME)/EpiFire/src
LDFLAGS = $(HOME)/EpiFire/src/*.o

all: epi_sim_abc

epi_sim_abc: epi_sim_abc.cpp ../MultiSeason_Network/MultiSeason_Sim.h
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) epi_sim_abc.cpp -o epi_sim_abc

epi_cluster_abc: epi_cluster_abc.cpp ../MultiSeason_Network/MultiSeason_Sim.h
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) epi_cluster_abc.cpp -o epi_cluster_abc

valgrind_test: valgrind_test.cpp ../MultiSeason_Network/MultiSeason_Sim.h
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) valgrind_test.cpp -o valgrind_test

test: test.cpp ../MultiSeason_Network/MultiSeason_Sim.h
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) test.cpp -o test 

