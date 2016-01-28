CC = g++
LINK = g++
AR = ar
CFLAGS = -O1 -g -Wall -fPIC
ARFLAGS = -rv
OUT_DIR = ./build
OBJS = $(OUT_DIR)/cheaptrick.o $(OUT_DIR)/common.o $(OUT_DIR)/d4c.o $(OUT_DIR)/dio.o $(OUT_DIR)/fft.o $(OUT_DIR)/matlabfunctions.o $(OUT_DIR)/stonemask.o $(OUT_DIR)/synthesis.o
LIBS =

default: $(OUT_DIR)/libworld.a
test: $(OUT_DIR)/test

$(OUT_DIR)/test: $(OUT_DIR)/libworld.a test/test.cpp test/audioio.cpp
	$(LINK) $(CFLAGS) -o $(OUT_DIR)/test test/test.cpp test/audioio.cpp $(OUT_DIR)/libworld.a -lm

$(OUT_DIR)/libworld.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libworld.a $(OBJS) $(LIBS)
	@echo Done.

$(OUT_DIR)/cheaptrick.o : src/cheaptrick.h src/common.h src/constantnumbers.h src/matlabfunctions.h
$(OUT_DIR)/common.o : src/common.h src/constantnumbers.h src/matlabfunctions.h
$(OUT_DIR)/d4c.o : src/d4c.h src/common.h src/constantnumbers.h src/matlabfunctions.h
$(OUT_DIR)/dio.o : src/dio.h src/common.h src/constantnumbers.h src/matlabfunctions.h
$(OUT_DIR)/fft.o : src/fft.h
$(OUT_DIR)/matlabfunctions.o : src/constantnumbers.h src/matlabfunctions.h
$(OUT_DIR)/stonemask.o : src/stonemask.h src/fft.h src/common.h src/constantnumbers.h src/matlabfunctions.h
$(OUT_DIR)/synthesis.o : src/synthesis.h src/common.h src/constantnumbers.h src/matlabfunctions.h

$(OUT_DIR)/%.o : src/%.cpp
	mkdir -p build
	$(CC) $(CFLAGS) -o $(OUT_DIR)/$*.o -c src/$*.cpp

clean:
	@echo 'Removing all temporary binaries... '
	@rm -f $(OUT_DIR)/libworld.a $(OUT_DIR)/*.o
	@echo Done.

clear:
	@echo 'Removing all temporary binaries... '
	@rm -f $(OUT_DIR)/libworld.a $(OUT_DIR)/*.o
	@echo Done.

