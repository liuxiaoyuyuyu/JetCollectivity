EXE=mymain06
GCC=g++

PYTHIAINC=${shell /home/fl52/Downloads/pythia8309/bin/pythia8-config --includedir}
PYTHIALIBS=${shell /home/fl52/Downloads/pythia8309/bin/pythia8-config --libdir}

FASTJETINC=${shell /home/fl52/Downloads/fastjet-install/bin/fastjet-config --cxxflags}
FASTJETLIBS=${shell /home/fl52/Downloads/fastjet-install/bin/fastjet-config --libs}

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

all: $(EXE)

mymain06.o: mymain06.cc
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) $(ROOTCFLAGS) -c mymain06.cc

$(EXE): mymain06.o
	$(GCC) -O3 -I$(PYTHIAINC) $(FASTJETINC)  $(ROOTCFLAGS) -o $(EXE) mymain06.o $(PYTHIALIBS)/libpythia8.a $(ROOTLIBS) -ldl $(FASTJETLIBS) -lRecursiveTools

.PHONY:clean
clean:
	-rm *.o $(EXE)
