CXX = g++
FOPENMP?=-fopenmp
CXXFLAGS = -Wall  -Wextra -std=c++0x $(FOPENMP) -I./gzstream -O2 -Wfatal-errors -D NDEBUG ${CPPFLAGS} 
#CXXFLAGS = -Wall -std=c++11 $(FOPENMP) -O3 -Wfatal-errors
PROGS1 = classify db_sort set_lcas db_shrink build_taxdb read_uid_mapping count_unique dump_taxdb 
TEST_PROGS = grade_classification test_hll_on_db dump_db_kmers
#PROGS = $(PROGS1) $(TEST_PROGS)
PROGS = $(PROGS1)
#LIBFLAGS = -L. -lz -lgzstream ${LDFLAGS}
LIBFLAGS = -L. -lz ${LDFLAGS}

.PHONY: all install clean

all: $(PROGS1)

allall: $(PROGS1) $(TEST_PROGS)

install: $(PROGS)
	cp $(PROGS) $(KRAKEN_DIR)/

clean:
	rm -rf $(PROGS) $(TEST_PROGS) *.o *.dSYM *.gch

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.o

grade_classification: taxdb.h report-cols.h

read_uid_mapping: quickfile.o

count_unique: hyperloglogplus.o

test_hll_on_db: krakendb.o hyperloglogplus.o quickfile.o

dump_db_kmers: krakendb.o quickfile.o

classify: classify.cpp krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.o gzstream.o hyperloglogplus.o taxdb.h report-cols.h readcounts.hpp
	$(CXX) $(CXXFLAGS) -o classify $^ $(LIBFLAGS)

build_taxdb: taxdb.h report-cols.h quickfile.o

make_seqid_to_taxid_map: quickfile.o

read_uid_mapping: quickfile.o krakenutil.o uid_mapping.o

krakenutil.o: krakenutil.cpp krakenutil.hpp taxdb.h report-cols.h
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c seqreader.cpp

gzstream.o: gzstream/gzstream.C gzstream/gzstream.h
	$(CXX) $(CXXFLAGS) -c -O gzstream/gzstream.C

quickfile.o: quickfile.cpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c quickfile.cpp

uid_mapping.o: krakenutil.hpp uid_mapping.hpp uid_mapping.cpp
	$(CXX) $(CXXFLAGS) -c uid_mapping.cpp

hyperloglogplus.o: hyperloglogplus.h hyperloglogplus.cpp
	$(CXX) $(CXXFLAGS) -c hyperloglogplus.cpp
