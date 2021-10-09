CXX = g++
FOPENMP?=-fopenmp
NDEBUG=-D NDEBUG

is_fopenmp_supported := $(shell touch foo.cpp && $(CXX) -fopenmp -c foo.cpp -o foo.o &> /dev/null && printf 'yes' || printf 'no'; rm -f foo.c foo.o)

$(info Compiling with multithreading support: $(is_fopenmp_supported))
ifeq ($(is_fopenmp_supported), no)
$(info To compile with OpenMP support, set CXX to a compiler that supports OpenMP, such as g++ or LLVM clang on OSX.)
FOPENMP=
endif


CXXFLAGS = -Wall -Wextra -Wfatal-errors -pipe -O2 -std=c++11 $(FOPENMP) -I./gzstream $(NDEBUG) ${CPPFLAGS}
PROGS1 = classify classifyExact db_sort set_lcas db_shrink build_taxdb read_uid_mapping count_unique dump_taxdb query_taxdb
TEST_PROGS = grade_classification test_hll_on_db dump_db_kmers
PROGS = $(PROGS1)
LDFLAGS = -L.
LDLIBS=-lz

.PHONY: all install clean debug

all: $(PROGS1)

debug: CXXFLAGS += -g -Og
debug: $(PROGS1)

allall: $(PROGS1) $(TEST_PROGS)

install: $(PROGS)
	install $(PROGS) $(KRAKEN_DIR)/

clean:
	rm -rf $(PROGS) $(TEST_PROGS) *.o *.dSYM *.gch

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: set_lcas.cpp krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.o gzstream.o

grade_classification: #taxdb.hpp report-cols.hpp

read_uid_mapping: quickfile.o

count_unique: hyperloglogplus.o seqreader.o krakenutil.o gzstream.o

test_count_unique: hyperloglogplus.o 

test_hll_on_db: krakendb.o hyperloglogplus.o quickfile.o

dump_db_kmers: krakendb.o quickfile.o

classify: classify.cpp krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.o gzstream.o hyperloglogplus.o

classifyExact: classify.cpp krakendb.o quickfile.o krakenutil.o seqreader.o uid_mapping.o gzstream.o hyperloglogplus.o
	$(CXX) $(CXXFLAGS) -DEXACT_COUNTING -o classifyExact $^ $(LDFLAGS) $(LDLIBS)

query_taxdb: #taxdb.hpp

build_taxdb: quickfile.o #taxdb.hpp report-cols.hpp

make_seqid_to_taxid_map: quickfile.o

read_uid_mapping: read_uid_mapping.cpp quickfile.o krakenutil.o uid_mapping.o gzstream.o

krakenutil.o: krakenutil.cpp krakenutil.hpp taxdb.hpp report-cols.hpp gzstream/gzstream.h

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp

gzstream.o: gzstream/gzstream.C gzstream/gzstream.h
	$(CXX) $(CXXFLAGS) -c -O gzstream/gzstream.C

quickfile.o: quickfile.cpp quickfile.hpp

uid_mapping.o: krakenutil.hpp uid_mapping.hpp uid_mapping.cpp

hyperloglogplus.o: hyperloglogplus.hpp hyperloglogplus.cpp
