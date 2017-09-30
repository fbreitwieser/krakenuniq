CXX = g++
FOPENMP?=-fopenmp
CXXFLAGS = -Wall -std=c++11 $(FOPENMP) -O2 -g -Wfatal-errors
PROGS = db_sort set_lcas classify make_seqid_to_taxid_map db_shrink build_taxdb grade_classification
LIBFLAGS = -L. -I./gzstream -L./gzstream -lz -lgzstream

.PHONY: all install clean

all: $(PROGS)

install: $(PROGS)
	cp $(PROGS) $(KRAKEN_DIR)/

clean:
	rm -f $(PROGS) *.o

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.cpp

grade_classification: taxdb.h

classify: krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.cpp
	$(CXX) $(CXXFLAGS) -o classify classify.cpp $^ $(LIBFLAGS)

build_taxdb: taxdb.h

make_seqid_to_taxid_map: quickfile.o

krakenutil.o: krakenutil.cpp krakenutil.hpp taxdb.h
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c seqreader.cpp

quickfile.o: quickfile.cpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c quickfile.cpp
