CXX = g++

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	CXXFLAGS = -std=c++0x -O3 -funroll-loops -Wall -Wno-sign-compare -fopenmp -finline-limit-50000 -lm -Dnullptr=0 -DSUB_BIN_PATH_DEFINED=\"$(PWD)/sub_bin\"
else
	CXXFLAGS = -std=c++0x -O3 -funroll-loops -Wno-sign-compare -Xpreprocessor -fopenmp -DSUB_BIN_PATH_DEFINED=\"$(PWD)/sub_bin\"
	LDFLAGS = -lomp
endif

PRG = platanus_b
OBJ = main.o assemble.o scaffold.o scaffoldGraph.o gapClose.o common.o baseCommand.o seqlib.o mapper.o gapCloseOLC.o merge.o iterate.o combine.o polish.o solveDBG.o pairedDBG.o phase.o consensus.o divide.o kmer_divide.o
SUBDIRS = minimap2
SUBBIN = sub_bin

.PHONY:all clean $(SUBDIRS) $(SUBBIN)


all: $(SUBDIRS) sub_bin $(PRG) 

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBBIN):
	mkdir -p sub_bin
	cp minimap2/minimap2 sub_bin
	cp scripts/*.pl scripts/*.pm  sub_bin

$(PRG): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
.cpp.o:
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -rf $(PRG) $(OBJ) $(SUBBIN)
	$(MAKE) clean -C $(SUBDIRS)
