CXX = g++
C99 = gcc -std=c99
LINK = g++
AR = ar
CXXFLAGS = -O1 -g -Wall -fPIC
CFLAGS = $(CXXFLAGS)
ARFLAGS = -rv
OUT_DIR = ./build
OBJS = $(OUT_DIR)/objs/cheaptrick.o $(OUT_DIR)/objs/common.o $(OUT_DIR)/objs/d4c.o $(OUT_DIR)/objs/dio.o $(OUT_DIR)/objs/fft.o $(OUT_DIR)/objs/matlabfunctions.o $(OUT_DIR)/objs/stonemask.o $(OUT_DIR)/objs/synthesis.o $(OUT_DIR)/objs/synthesisrealtime.o
LIBS =

default: $(OUT_DIR)/libworld.a
test: $(OUT_DIR)/test $(OUT_DIR)/ctest

test_OBJS=$(OUT_DIR)/objs/test/audioio.o $(OUT_DIR)/objs/test/test.o
$(OUT_DIR)/test: $(OUT_DIR)/libworld.a $(test_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/test $(test_OBJS) $(OUT_DIR)/libworld.a -lm

ctest_OBJS=$(OUT_DIR)/objs/test/audioio.o $(OUT_DIR)/objs/test/ctest.o
$(OUT_DIR)/ctest: $(OUT_DIR)/libworld.a $(ctest_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/ctest $(ctest_OBJS) $(OUT_DIR)/libworld.a -lm

$(OUT_DIR)/libworld.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libworld.a $(OBJS) $(LIBS)
	@echo Done.

$(OUT_DIR)/objs/test/audioio.o : test/audioio.h
$(OUT_DIR)/objs/test/test.o : test/audioio.h src/world/d4c.h src/world/dio.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/test/ctest.o : test/audioio.h src/world/d4c.h src/world/dio.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/cheaptrick.o : src/world/cheaptrick.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/common.o : src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/d4c.o : src/world/d4c.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/dio.o : src/world/dio.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/fft.o : src/world/fft.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/matlabfunctions.o : src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/stonemask.o : src/world/stonemask.h src/world/fft.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/synthesis.o : src/world/synthesis.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/synthesisrealtime.o : src/world/synthesisrealtime.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h

$(OUT_DIR)/objs/test/%.o : test/%.c
	mkdir -p $(OUT_DIR)/objs/test
	$(C99) $(CFLAGS) -Isrc -o "$@" -c "$<"

$(OUT_DIR)/objs/%.o : src/%.c
	mkdir -p $(OUT_DIR)/objs
	$(C99) $(CFLAGS) -Isrc -o "$@" -c "$<"

$(OUT_DIR)/objs/test/%.o : test/%.cpp
	mkdir -p $(OUT_DIR)/objs/test
	$(CXX) $(CXXFLAGS) -Isrc -o "$@" -c "$<"

$(OUT_DIR)/objs/%.o : src/%.cpp
	mkdir -p $(OUT_DIR)/objs
	$(CXX) $(CXXFLAGS) -Isrc -o "$@" -c "$<"

clean:
	@echo 'Removing all temporary binaries... '
	@$(RM) $(OUT_DIR)/libworld.a $(OBJS)
	@$(RM) $(test_OBJS) $(ctest_OBJS) $(OUT_DIR)/test $(OUT_DIR)/ctest
	@echo Done.

clear: clean

.PHONY: clean clear test default
.DELETE_ON_ERRORS: