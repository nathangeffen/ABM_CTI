BASE=abm_spaces
EXE=$(BASE)
CORES=0
ifeq ($(OS),Windows_NT)
	EXE=$(BASE)+.exe
endif

release:
	g++ -O3 -std=c++11 -Wall -pedantic -DNDEBUG $(BASE).cc -o $(EXE) -lpthread

debug:
	g++ -g -std=c++11 -Wall -pedantic abm_spaces.cc -o $(EXE) -lpthread

analyze10000five: release
	./$(EXE) --num_runs=10000 --cores=$(CORES) >a.csv
	./$(EXE) --num_runs=10000 --tracing=1 --cores=$(CORES) > b.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=8 --cores=$(CORES) > c.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --cores=$(CORES) > d.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=8 --cores=$(CORES) > e.csv
	R -f analyze5.R

analyzeTaT2: release
	./$(EXE) --num_runs=10000 --tracing=1 --tat=2 --cores=$(CORES) > t2.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=3 --cores=$(CORES) > t3.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=4 --cores=$(CORES) > t4.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=5 --cores=$(CORES) > t5.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=6 --cores=$(CORES) > t6.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=7 --cores=$(CORES) > t7.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=8 --cores=$(CORES) > t8.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=9 --cores=$(CORES) > t9.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=10 --cores=$(CORES) > t10.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=11 --cores=$(CORES) > t11.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=2 --cores=$(CORES) > t2n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=3 --cores=$(CORES) > t3n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=4 --cores=$(CORES) > t4n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=5 --cores=$(CORES) > t5n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=6 --cores=$(CORES) > t6n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=7 --cores=$(CORES) > t7n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=8 --cores=$(CORES) > t8n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=9 --cores=$(CORES) > t9n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=10 --cores=$(CORES) > t10n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=11 --cores=$(CORES) > t11n.csv
	R -f analyzeTaT2.R

analyzeUltra: release
	./$(EXE) --num_runs=10000 --tracing=0 --tat=2 --cores=$(CORES) > t2x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=3 --cores=$(CORES) > t3x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=4 --cores=$(CORES) > t4x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=5 --cores=$(CORES) > t5x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=6 --cores=$(CORES) > t6x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=7 --cores=$(CORES) > t7x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=8 --cores=$(CORES) > t8x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=9 --cores=$(CORES) > t9x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=10 --cores=$(CORES) > t10x.csv
	./$(EXE) --num_runs=10000 --tracing=0 --tat=11 --cores=$(CORES) > t11x.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=2 --cores=$(CORES) > t2.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=3 --cores=$(CORES) > t3.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=4 --cores=$(CORES) > t4.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=5 --cores=$(CORES) > t5.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=6 --cores=$(CORES) > t6.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=7 --cores=$(CORES) > t7.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=8 --cores=$(CORES) > t8.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=9 --cores=$(CORES) > t9.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=10 --cores=$(CORES) > t10.csv
	./$(EXE) --num_runs=10000 --tracing=1 --tat=11 --cores=$(CORES) > t11.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=2 --cores=$(CORES) > t2n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=3 --cores=$(CORES) > t3n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=4 --cores=$(CORES) > t4n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=5 --cores=$(CORES) > t5n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=6 --cores=$(CORES) > t6n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=7 --cores=$(CORES) > t7n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=8 --cores=$(CORES) > t8n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=9 --cores=$(CORES) > t9n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=10 --cores=$(CORES) > t10n.csv
	./$(EXE) --num_runs=10000 --test_contacts=0 --tracing=1 --tat=11 --cores=$(CORES) > t11n.csv
	R -f analyzeUltra.R


abm_hetgen: abm_hetgen.cc
	g++ -std=c++17 -O3 -DNDEBUG -Wall -o abm_hetgen abm_hetgen.cc -lpthread

abm_hetgen_debug: abm_hetgen.cc
	g++ -std=c++17 -g -Wall -o abm_hetgen_debug abm_hetgen.cc -lpthread
