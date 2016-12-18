CXX = g++
C99 = gcc -std=c99
LINK = g++
AR = ar
#DEBUG_FLAG=-g
CXXFLAGS = -O1 -Wall -fPIC $(DEBUG_FLAG)
CFLAGS = $(CXXFLAGS)
ARFLAGS = -rv
OUT_DIR = ./build
OBJS = $(OUT_DIR)/objs/world/cheaptrick.o $(OUT_DIR)/objs/world/common.o $(OUT_DIR)/objs/world/d4c.o $(OUT_DIR)/objs/world/dio.o $(OUT_DIR)/objs/world/fft.o $(OUT_DIR)/objs/world/harvest.o $(OUT_DIR)/objs/world/matlabfunctions.o $(OUT_DIR)/objs/world/stonemask.o $(OUT_DIR)/objs/world/synthesis.o $(OUT_DIR)/objs/world/synthesisrealtime.o
LIBS =

all: default test

###############################################################################################################
### Tests
###############################################################################################################
test: $(OUT_DIR)/test $(OUT_DIR)/ctest

test_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/test.o
$(OUT_DIR)/test: $(OUT_DIR)/libworld.a $(test_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/test $(test_OBJS) $(OUT_DIR)/libworld.a -lm

ctest_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/ctest.o
$(OUT_DIR)/ctest: $(OUT_DIR)/libworld.a $(ctest_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/ctest $(ctest_OBJS) $(OUT_DIR)/libworld.a -lm

$(OUT_DIR)/objs/test/test.o : tools/audioio.h src/world/d4c.h src/world/dio.h src/world/harvest.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/test/ctest.o : tools/audioio.h src/world/d4c.h src/world/dio.h src/world/harvest.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h

###############################################################################################################
### Library
###############################################################################################################
default: $(OUT_DIR)/libworld.a

$(OUT_DIR)/libworld.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libworld.a $(OBJS) $(LIBS)
	@echo Done.

$(OUT_DIR)/objs/world/cheaptrick.o : src/world/cheaptrick.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/common.o : src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/d4c.o : src/world/d4c.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/dio.o : src/world/dio.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/fft.o : src/world/fft.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/harvest.o : src/world/harvest.h src/world/fft.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/matlabfunctions.o : src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/stonemask.o : src/world/stonemask.h src/world/fft.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/synthesis.o : src/world/synthesis.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/synthesisrealtime.o : src/world/synthesisrealtime.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h


###############################################################################################################
### Global rules
###############################################################################################################
$(OUT_DIR)/objs/test/%.o : test/%.c
	mkdir -p $(OUT_DIR)/objs/test
	$(C99) $(CFLAGS) -Isrc -Itools -o "$@" -c "$<"

$(OUT_DIR)/objs/test/%.o : test/%.cpp
	mkdir -p $(OUT_DIR)/objs/test
	$(CXX) $(CXXFLAGS) -Isrc -Itools -o "$@" -c "$<"

$(OUT_DIR)/objs/tools/%.o : tools/%.cpp
	mkdir -p $(OUT_DIR)/objs/tools
	$(CXX) $(CXXFLAGS) -Isrc -o "$@" -c "$<"

$(OUT_DIR)/objs/world/%.o : src/%.cpp
	mkdir -p $(OUT_DIR)/objs/world
	$(CXX) $(CXXFLAGS) -Isrc -o "$@" -c "$<"

clean:
	@echo 'Removing all temporary binaries... '
	@$(RM) $(OUT_DIR)/libworld.a $(OBJS)
	@$(RM) $(test_OBJS) $(ctest_OBJS) $(OUT_DIR)/test $(OUT_DIR)/ctest
	@echo Done.

clear: clean

.PHONY: clean clear test default
.DELETE_ON_ERRORS:
