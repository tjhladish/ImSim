INCLUDE = -I /home/tjhladish/work/EpiFire/src
LDFLAGS = /home/tjhladish/work/EpiFire/src/*.o

all: tmp_research

tmp_research: tmp_research.cpp
	g++ $(INCLUDE) $(LDFLAGS) tmp_research.cpp -o tmp_research

epi_sim: epi_sim.cpp MultiSeason_Sim.h
	g++ -O2 $(INCLUDE) $(LDFLAGS) MultiSeason_Sim.h epi_sim.cpp -o epi_sim

