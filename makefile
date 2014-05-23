# ----------------------------------------------------
# CPLEX Solver makefile
# Pedro Dias
# v0.1
# ----------------------------------------------------

# System Information
SYSTEM    = x86-64_sles10_4.1
LIBFORMAT = static_pic

# User Information
USER = $(shell whoami)
HOME = /home/$(USER)

# CPLEX Directories Information
CPLEXDIR    = $(HOME)/CPLEX/cplex
CPLEXBINDIR = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXINCDIR = $(CPLEXDIR)/include

# Compiler Selections and Options
CC   = gcc -O0
COPT = -m64 -fPIC -fno-strict-aliasing

# Solver Directories
SRC    = ./src
INC    = ./include
BUILD  = ./build
TARGET = ./bin/solver

# Make Options
SRCEXT  = c
SOURCES = $(shell find $(SRC) -type f -name *.$(SRCEXT))
OBJECTS = $(patsubst $(SRC)/%, $(BUILD)/%, $(SOURCES:.$(SRCEXT)=.o))

# Link Options
CLNDIRS  = -L$(CPLEXLIBDIR)
CLNFLAGS = -lcplex -lm -pthread
CFLAGS   = $(COPT) -I$(CPLEXINCDIR) -I$(INC)

###############
#### Rules ####
###############

$(TARGET): $(OBJECTS)
	@echo "Linking... "
	$(CC) $(CFLAGS) $(CLNDIRS) $^ -o $(TARGET) $(CLNFLAGS)

$(BUILD)/main.o: $(SRC)/main.c
	@echo "Compiling src/main.c... "
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD)/network.o: $(SRC)/network.c
	@echo "Compiling src/network.c... "
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	@echo "Cleaning... "
	rm -vf ./build/*.o
	rm -vf ./bin/*
