CXX = g++
FOPENMP?=-fopenmp
CXXFLAGS = -Wall -std=c++11 $(FOPENMP) -O2 -g -Wfatal-errors
#CXXFLAGS = -Wall -std=c++11 $(FOPENMP) -O3 -Wfatal-errors
PROGS = db_sort set_lcas classify make_seqid_to_taxid_map db_shrink build_taxdb grade_classification dump_taxdb read_uid_mapping
LIBFLAGS = -L. -I./gzstream -L./gzstream -lz -lgzstream

.PHONY: all install clean

all: $(PROGS)

install: $(PROGS)
	cp $(PROGS) $(KRAKEN_DIR)/

clean:
	rm -f $(PROGS) *.o

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.o

grade_classification: taxdb.h

read_uid_mapping: quickfile.o

classify: classify.cpp krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.o hyperloglogplus.h
	$(CXX) $(CXXFLAGS) -o classify $^ $(LIBFLAGS)

build_taxdb: taxdb.h

make_seqid_to_taxid_map: quickfile.o

read_uid_mapping: quickfile.o krakenutil.o uid_mapping.o

krakenutil.o: krakenutil.cpp krakenutil.hpp taxdb.h
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c seqreader.cpp

quickfile.o: quickfile.cpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c quickfile.cpp

uid_mapping.o: krakenutil.hpp uid_mapping.hpp uid_mapping.cpp
	$(CXX) $(CXXFLAGS) -c uid_mapping.cpp
