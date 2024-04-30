SRC_DIR=src

INCLUDE_DIR=include

CXX=g++

ifdef D
	OPT=-g -O0
else
	OPT=-O3 -DNDEBUG
endif

CXXFLAGS=$(OPT) -I$(INCLUDE_DIR)
LDFLAGS=-lssl -lcrypto

ALL=unif_gen adv_gen

all: $(ALL)

unif_gen: $(SRC_DIR)/unif_gen.cpp
	$(CXX) $(CXXFLAGS) -o unif_gen $(SRC_DIR)/unif_gen.cpp $(LDFLAGS)

adv_gen: $(SRC_DIR)/adv_gen.cpp
	$(CXX) $(CXXFLAGS) -o adv_gen $(SRC_DIR)/adv_gen.cpp $(LDFLAGS)

clean:
	rm -f $(ALL)
