ROOT_DIR=$(shell pwd)
ODIR  = $(ROOT_DIR)/obj
SDIR  = $(ROOT_DIR)/src

CXX   = g++
CFLAG = -std=c++11 -lfftw3
DFLAG = -D__SPECTRA__
 
DEPS  = $(shell ls $(SDIR)/*.h)
SRC   = $(shell ls $(SDIR)/*.cpp)
OBJ   = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRC))

analyser.x : $(OBJ)
	$(CXX) -o $@ $^ $(DFLAG) $(CFLAG)

$(ODIR)/%.o : $(SDIR)/%.cpp $(DEPS) | $(ODIR)/.
	$(CXX) -c -o $@ $< $(CFLAG)

%/. : 
	mkdir -p $(patsubst %/.,%,$@)
	
.PRECIOUS: %/.
.PHONY: clean clean_all
clean:
	rm -rf analyser.x $(ODIR)
clean_all:
	rm -rf analyser.x *.dat $(ODIR)
