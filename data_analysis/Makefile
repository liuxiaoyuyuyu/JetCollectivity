EXE=new_default_data_vn
GCC=g++

PYTHIAINC=${shell /afs/cern.ch/user/x/xiaoyul/pythia8312/bin/pythia8-config --includedir}
PYTHIALIBS=${shell /afs/cern.ch/user/x/xiaoyul/pythia8312/bin/pythia8-config --libdir}

FASTJETINC=${shell /afs/cern.ch/user/x/xiaoyul/fastjet-install/bin/fastjet-config --cxxflags}
FASTJETLIBS=${shell /afs/cern.ch/user/x/xiaoyul/fastjet-install/bin/fastjet-config --libs}

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

all: $(EXE)

new_default_data_vn.o: new_default_data_vn.cc
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) $(ROOTCFLAGS) -c new_default_data_vn.cc

$(EXE): new_default_data_vn.o
	$(GCC) -O3 -I$(PYTHIAINC) $(FASTJETINC)  $(ROOTCFLAGS) -o $(EXE) new_default_data_vn.o $(PYTHIALIBS)/libpythia8.a $(ROOTLIBS) -ldl $(FASTJETLIBS) -lRecursiveTools

.PHONY:clean
clean:
	-rm *.o $(EXE)
