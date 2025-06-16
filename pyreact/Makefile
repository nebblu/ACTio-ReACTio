COPTER_DIR = ../reactions
COPTER_LIB = $(COPTER_DIR)/lib/libcopter.a

WRAPPER_LIB = libreact_wrapper.so

# Compiler settings
CC=/opt/homebrew/opt/llvm/bin/clang
CXX=/opt/homebrew/opt/llvm/bin/clang++

CFLAGS += -fPIC -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
CXXFLAGS += -Wall -fPIC -std=c++11 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include

# Linking flags
LDFLAGS += -L/opt/homebrew/opt/gsl/lib \
           -L/Users/bbose/Desktop/cosmo_codes/sundials-4.1.0/instdir/lib \
           -L/opt/homebrew/opt/libomp/lib \
           -L/opt/homebrew/opt/llvm/lib \
           -Wl,-rpath,/Users/bbose/Desktop/cosmo_codes/sundials-4.1.0/instdir/lib

LDFLAGS += -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -lomp

CPPFLAGS += -I/opt/homebrew/opt/gsl/include \
            -I/Users/bbose/Desktop/cosmo_codes/sundials-4.1.0/instdir/include \
            -I/opt/homebrew/opt/llvm/include

all: $(WRAPPER_LIB)

$(WRAPPER_LIB): react_wrapper.cpp $(COPTER_LIB)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -shared -I$(COPTER_DIR)/include

$(COPTER_LIB):
	@echo "\nBuilding COPTER\n"
	cd $(COPTER_DIR) && CXX=$(CXX) CC=$(CC) CXXFLAGS="$(CXXFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)" ./build_copter.sh

# Clean up
.PHONY: clean
clean:
	rm -f $(WRAPPER_LIB)
	rm -rf $(WRAPPER_LIB).dSYM/
	@echo "Entering $(COPTER_DIR) and running make clean.\n"
	cd $(COPTER_DIR) && make clean && make uninstall && rm Makefile