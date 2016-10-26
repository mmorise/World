CXX = g++
C99 = gcc -std=c99
LINK = g++
AR = ar
#DEBUG_FLAG=-g
CXXFLAGS = -O1 -Wall -fPIC $(DEBUG_FLAG)
CFLAGS = $(CXXFLAGS)
ARFLAGS = -rv
OUT_DIR = ./../build
OBJS = $(OUT_DIR)/objs/world/cheaptrick.o $(OUT_DIR)/objs/world/common.o $(OUT_DIR)/objs/world/d4c.o $(OUT_DIR)/objs/world/dio.o $(OUT_DIR)/objs/world/fft.o $(OUT_DIR)/objs/world/matlabfunctions.o $(OUT_DIR)/objs/world/stonemask.o $(OUT_DIR)/objs/world/synthesis.o $(OUT_DIR)/objs/world/synthesisrealtime.o
LIBS =

all: default example

default: $(OUT_DIR)/libworld.a
###############################################################################################################
### Examples
### If you want to add your example, please add the compile command.
### Note: You should write the delete command in "clean:" section.
###############################################################################################################
examples: $(OUT_DIR)/analysis $(OUT_DIR)/synthesis

$(OUT_DIR)/analysis: $(OUT_DIR)/objs/examples/analysis.o $(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/libworld.a
	 $(LINK) $(CXXFLAGS) -o $@ $^ -lm

$(OUT_DIR)/synthesis: $(OUT_DIR)/objs/examples/synthesis.o $(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/libworld.a
	$(LINK) $(CXXFLAGS) -o $@ $^ -lm

$(OUT_DIR)/objs/tools/audioio.o : ../tools/audioio.h

###############################################################################################################
### Library
###############################################################################################################
$(OUT_DIR)/libworld.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libworld.a $(OBJS) $(LIBS)
	@echo Done.

$(OUT_DIR)/objs/world/cheaptrick.o : ../src/world/cheaptrick.h ../src/world/common.h ../src/world/constantnumbers.h ../src/world/matlabfunctions.h ../src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/common.o : ../src/world/common.h ../src/world/constantnumbers.h ../src/world/matlabfunctions.h ../src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/d4c.o : ../src/world/d4c.h ../src/world/common.h ../src/world/constantnumbers.h ../src/world/matlabfunctions.h ../src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/dio.o : ../src/world/dio.h ../src/world/common.h ../src/world/constantnumbers.h ../src/world/matlabfunctions.h ../src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/fft.o : ../src/world/fft.h ../src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/matlabfunctions.o : ../src/world/constantnumbers.h ../src/world/matlabfunctions.h ../src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/stonemask.o : ../src/world/stonemask.h ../src/world/fft.h ../src/world/common.h ../src/world/constantnumbers.h ../src/world/matlabfunctions.h ../src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/synthesis.o : ../src/world/synthesis.h ../src/world/common.h ../src/world/constantnumbers.h ../src/world/matlabfunctions.h ../src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/synthesisrealtime.o : ../src/world/synthesisrealtime.h ../src/world/common.h ../src/world/constantnumbers.h ../src/world/matlabfunctions.h ../src/world/macrodefinitions.h


###############################################################################################################
### Global rules
###############################################################################################################
$(OUT_DIR)/objs/examples/%.o : ../examples/%.cpp
	mkdir -p $(OUT_DIR)/objs/examples
	$(CXX) $(CXXFLAGS) -I../src -I../tools  -o "$@" -c "$<"

$(OUT_DIR)/objs/tools/%.o : ../tools/%.cpp
	mkdir -p $(OUT_DIR)/objs/tools
	$(CXX) $(CXXFLAGS) -I../src -o "$@" -c "$<"

$(OUT_DIR)/objs/world/%.o : ../src/%.cpp
	mkdir -p $(OUT_DIR)/objs/world
	$(CXX) $(CXXFLAGS) -I../src -o "$@" -c "$<"

clean:
	@echo 'Removing all temporary binaries... '
	@$(RM) $(OUT_DIR)/libworld.a $(OBJS)
	@$(RM) $(OUT_DIR)/analysis $(OUT_DIR)/synthesis
	@echo Done.

clear: clean

.PHONY: clean clear test default
.DELETE_ON_ERRORS:
