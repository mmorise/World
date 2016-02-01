CXX = g++
C99 = gcc -std=c99
LINK = g++
AR = ar
CXXFLAGS = -O1 -g -Wall -fPIC
CFLAGS = $(CXXFLAGS)
ARFLAGS = -rv
OUT_DIR = ./build
OBJS = $(OUT_DIR)/objs/cheaptrick.o $(OUT_DIR)/objs/common.o $(OUT_DIR)/objs/d4c.o $(OUT_DIR)/objs/dio.o $(OUT_DIR)/objs/fft.o $(OUT_DIR)/objs/matlabfunctions.o $(OUT_DIR)/objs/stonemask.o $(OUT_DIR)/objs/synthesis.o
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
$(OUT_DIR)/objs/test/test.o : test/audioio.h src/d4c.h src/dio.h src/matlabfunctions.h src/cheaptrick.h src/stonemask.h src/synthesis.h src/common.h src/fft.h src/macrodefinitions.h
$(OUT_DIR)/objs/test/ctest.o : test/audioio.h src/d4c.h src/dio.h src/matlabfunctions.h src/cheaptrick.h src/stonemask.h src/synthesis.h src/common.h src/fft.h src/macrodefinitions.h
$(OUT_DIR)/objs/cheaptrick.o : src/cheaptrick.h src/common.h src/constantnumbers.h src/matlabfunctions.h src/macrodefinitions.h
$(OUT_DIR)/objs/common.o : src/common.h src/constantnumbers.h src/matlabfunctions.h src/macrodefinitions.h
$(OUT_DIR)/objs/d4c.o : src/d4c.h src/common.h src/constantnumbers.h src/matlabfunctions.h src/macrodefinitions.h
$(OUT_DIR)/objs/dio.o : src/dio.h src/common.h src/constantnumbers.h src/matlabfunctions.h src/macrodefinitions.h
$(OUT_DIR)/objs/fft.o : src/fft.h src/macrodefinitions.h
$(OUT_DIR)/objs/matlabfunctions.o : src/constantnumbers.h src/matlabfunctions.h src/macrodefinitions.h
$(OUT_DIR)/objs/stonemask.o : src/stonemask.h src/fft.h src/common.h src/constantnumbers.h src/matlabfunctions.h src/macrodefinitions.h
$(OUT_DIR)/objs/synthesis.o : src/synthesis.h src/common.h src/constantnumbers.h src/matlabfunctions.h src/macrodefinitions.h

$(OUT_DIR)/objs/test/%.o : test/%.c
	mkdir -p $(OUT_DIR)/objs/test
	$(C99) $(CFLAGS) -o "$@" -c "$<"

$(OUT_DIR)/objs/%.o : src/%.c
	mkdir -p $(OUT_DIR)/objs
	$(C99) $(CFLAGS) -o "$@" -c "$<"

$(OUT_DIR)/objs/test/%.o : test/%.cpp
	mkdir -p $(OUT_DIR)/objs/test
	$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OUT_DIR)/objs/%.o : src/%.cpp
	mkdir -p $(OUT_DIR)/objs
	$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

clean:
	@echo 'Removing all temporary binaries... '
	@$(RM) $(OUT_DIR)/libworld.a $(OBJS)
	@$(RM) $(test_OBJS) $(ctest_OBJS) $(OUT_DIR)/test $(OUT_DIR)/ctest
	@echo Done.

clear: clean

.PHONY: clean clear
.DELETE_ON_ERRORS:
