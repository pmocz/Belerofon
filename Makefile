EXEC        = belerofon
CONFIG      = Config.sh
BUILD_DIR   = build
SRC_DIR     = src

ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
else
SYSTYPE := "default"
endif

$(info Build configuration:)
$(info SYSTYPE: $(SYSTYPE))



CC          = gcc
CPP         = gcc
OPTIMIZE    = -std=c99 -fopenmp -O3 -g -Wall -Wno-unknown-pragmas
HDF5_INCL   = -I/usr/include/hdf5/serial
HDF5_LIB    = -L/usr/lib/hdf5/serial  -lhdf5_serial -lsz -lz -ldl -lm 
FFTW_INCL   = -I/usr/local/include
FFTW_LIB    = -L/usr/local/lib  -lfftw3_omp -lfftw3 -lm
LINKER      = gcc


ifeq ($(SYSTYPE),"tiger")
CC          = mpicc
CPP         = mpicc
OPTIMIZE    = -std=c99 -fopenmp -O3 -m64 -g -Wall -Wno-unknown-pragmas
HDF5_INCL   = 
HDF5_LIB    = -L$(LD_LIBRARY_PATH) -lhdf5 -lsz -lz -ldl -lm
FFTW_INCL   = 
FFTW_LIB    = -lfftw3_omp -lfftw3 -lm
LINKER      = mpicc
endif


############################################
# needed objects/headers

SUBDIRS = .
OBJS = main.o utils.o
INCL = utils.h allvars.h



############################################
# combine options

CFLAGS = $(OPTIMIZE) $(HDF5_INCL) $(FFTW_INCL)

LIBS = $(HDF5_LIB) $(FFTW_LIB)




SUBDIRS := $(addprefix $(BUILD_DIR)/,$(SUBDIRS))
OBJS := $(addprefix $(BUILD_DIR)/,$(OBJS))
INCL := $(addprefix $(SRC_DIR)/,$(INCL))



############################################
# create subdirs

RESULT := $(shell mkdir -p $(SUBDIRS)  )



############################################
# build rules

all: $(EXEC)

$(EXEC): $(OBJS)
	$(LINKER) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -f $(OBJS) $(EXEC)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCL) $(MAKEFILES)
	$(CPP) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCL) $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

