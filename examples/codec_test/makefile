CXX = g++
C99 = gcc -std=c99
LINK = g++
AR = ar
#DEBUG_FLAG=-g
CXXFLAGS = -O1 -Wall -fPIC $(DEBUG_FLAG)
CFLAGS = $(CXXFLAGS)
ARFLAGS = -rv
OUT_DIR = ./../../build
IN_DIR = ./../..
OBJS = $(OUT_DIR)/objs/world/cheaptrick.o $(OUT_DIR)/objs/world/common.o $(OUT_DIR)/objs/world/d4c.o $(OUT_DIR)/objs/world/dio.o $(OUT_DIR)/objs/world/fft.o $(OUT_DIR)/objs/world/matlabfunctions.o $(OUT_DIR)/objs/world/stonemask.o $(OUT_DIR)/objs/world/synthesis.o $(OUT_DIR)/objs/world/synthesisrealtime.o $(OUT_DIR)/objs/world/harvest.o $(OUT_DIR)/objs/world/codec.o
LIBS =

all: default examples

default: $(OUT_DIR)/libworld.a
###############################################################################################################
### Examples
### If you want to add your example, please add the compile command.
### Note: You should write the delete command in "clean:" section.
###############################################################################################################
examples: $(OUT_DIR)/f0analysis $(OUT_DIR)/spanalysis $(OUT_DIR)/apanalysis $(OUT_DIR)/readandsynthesis

$(OUT_DIR)/f0analysis: $(OUT_DIR)/objs/examples/f0analysis.o $(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/tools/parameterio.o $(OUT_DIR)/libworld.a
	 $(LINK) $(CXXFLAGS) -o $@ $^ -lm

$(OUT_DIR)/spanalysis: $(OUT_DIR)/objs/examples/spanalysis.o $(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/tools/parameterio.o $(OUT_DIR)/libworld.a
	 $(LINK) $(CXXFLAGS) -o $@ $^ -lm

$(OUT_DIR)/apanalysis: $(OUT_DIR)/objs/examples/apanalysis.o $(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/tools/parameterio.o $(OUT_DIR)/libworld.a
	 $(LINK) $(CXXFLAGS) -o $@ $^ -lm

$(OUT_DIR)/readandsynthesis: $(OUT_DIR)/objs/examples/readandsynthesis.o $(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/tools/parameterio.o $(OUT_DIR)/libworld.a
	 $(LINK) $(CXXFLAGS) -o $@ $^ -lm

$(OUT_DIR)/objs/tools/audioio.o : $(IN_DIR)/tools/audioio.h
$(OUT_DIR)/objs/tools/parameterio.o : $(IN_DIR)/tools/parameterio.h

###############################################################################################################
### Library
###############################################################################################################
$(OUT_DIR)/libworld.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libworld.a $(OBJS) $(LIBS)
	@echo Done.

$(OUT_DIR)/objs/world/cheaptrick.o : $(IN_DIR)/src/world/cheaptrick.h $(IN_DIR)/src/world/common.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/common.o : $(IN_DIR)/src/world/common.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/d4c.o : $(IN_DIR)/src/world/d4c.h $(IN_DIR)/src/world/common.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/dio.o : $(IN_DIR)/src/world/dio.h $(IN_DIR)/src/world/common.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/fft.o : $(IN_DIR)/src/world/fft.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/matlabfunctions.o : $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/stonemask.o : $(IN_DIR)/src/world/stonemask.h $(IN_DIR)/src/world/fft.h $(IN_DIR)/src/world/common.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/synthesis.o : $(IN_DIR)/src/world/synthesis.h $(IN_DIR)/src/world/common.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/synthesisrealtime.o : $(IN_DIR)/src/world/synthesisrealtime.h $(IN_DIR)/src/world/common.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/harvest.o : $(IN_DIR)/src/world/harvest.h $(IN_DIR)/src/world/common.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/codec.o : $(IN_DIR)/src/world/codec.h $(IN_DIR)/src/world/constantnumbers.h $(IN_DIR)/src/world/matlabfunctions.h $(IN_DIR)/src/world/macrodefinitions.h

###############################################################################################################
### Global rules
###############################################################################################################
$(OUT_DIR)/objs/examples/%.o : $(IN_DIR)/examples/codec_test/%.cpp
	mkdir -p $(OUT_DIR)/objs/examples
	$(CXX) $(CXXFLAGS) -I$(IN_DIR)/src -I$(IN_DIR)/tools  -o "$@" -c "$<"

$(OUT_DIR)/objs/tools/%.o : $(IN_DIR)/tools/%.cpp
	mkdir -p $(OUT_DIR)/objs/tools
	$(CXX) $(CXXFLAGS) -I$(IN_DIR)/src -o "$@" -c "$<"

$(OUT_DIR)/objs/world/%.o : $(IN_DIR)/src/%.cpp
	mkdir -p $(OUT_DIR)/objs/world
	$(CXX) $(CXXFLAGS) -I$(IN_DIR)/src -o "$@" -c "$<"

clean:
	@echo 'Removing all temporary binaries... '
	@$(RM) $(OUT_DIR)/libworld.a $(OBJS)
	@$(RM) $(OUT_DIR)/f0analysis $(OUT_DIR)/spanalysis $(OUT_DIR)/apanalysis $(OUT_DIR)/readandsynthesis
	@echo Done.

clear: clean

.PHONY: clean clear test default
.DELETE_ON_ERRORS:
