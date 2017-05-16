#
# Makefile
#
# A Series of Randmoness Tests for Binary Sequence Validation
# by snovvcrash
# 04.2017
#

CXXTARGET = randomness_test
CTARGET   =

CXX = g++
CC  = gcc

CXXFLAGS += -Wall -c -std=c++11 -O2
CFLAGS   += -Wall -c
LDFLAGS  += -Wall

GSL_INCLUDE +=
GSL_LIBRARY +=

HEADERS    = $(wildcard *.hxx)
CXXSOURCES = $(wildcard *.cxx)
CXXOBJECTS = $(patsubst %.cxx, %.o, $(CXXSOURCES))
CSOURCES   = $(wildcard *.c)
COBJECTS   = $(patsubst %.c, %.o, $(CSOURCES))

.PHONY: cxxbuild cbuild all default clean
.PRECIOUS: $(CXXTARGET) $(CTARGET) $(CXXOBJECTS) $(COBJECTS)

all: clean default
default: cxxbuild
cxxbuild: $(CXXTARGET)
	@echo "Build cxx-project"
cbuild: $(CTARGET)
	@echo "Build c-project"

$(CXXTARGET): $(CXXOBJECTS)
	@echo "(CXX) $?"
	@$(CXX) $(CXXOBJECTS) $(LDFLAGS) $(GSL_LIBRARY) -o $@

$(CTARGET): $(COBJECTS)
	@echo "(CC) $?"
	@$(CC) $(COBJECTS) $(LDFLAGS) $(GSL_LIBRARY) -o $@
	
%.o: %.cxx $(HEADERS)
	@echo "(CXX) $<"
	@$(CXX) $(CXXFLAGS) $(GSL_INCLUDE) $< -o $@

%.o: %.c $(HEADERS)
	@echo "(CC) $<"
	@$(CC) $(CFLAGS) $(GSL_INCLUDE) $< -o $@

debug: CXXFLAGS += -DDEBUG -g -O0
debug: CFLAGS   += -DDEBUG -g -O0
debug: all
	@echo "DEBUG MODE"

clean:
	@echo "Clean project"
	@rm -rfv *.o $(CXXTARGET) $(CTARGET)
