CXX=g++
CCFLAGS= -g -std=c++14 -pthread
LDFLAGS= -lntl -lgmpxx -lgmp -lm
WARNFLAGS= -Wall -Wunknown-pragmas
OPTIMIZERS= -march=native -O3 -fomit-frame-pointer

MACHINE_PREFIX=/path/to/project

INCLUDE_DIRS= -I./ -I./.. -I./../dgs
LDFLAGS_DGS= -L$(MACHINE_PREFIX)/src/dgs/.libs -Wl,-rpath,$(MACHINE_PREFIX)/src/dgs/.libs -ldgs

BINDIR=bin
TARGET=$(BINDIR)/lwe-ot

SRCS= $(wildcard */*.cpp) $(wildcard *.cpp)
ASMS= $(wildcard prng/*.s)
OBJS= $(SRCS:.cpp=.o) $(ASMS:.s=.o)
HDRS= $(wildcard *.h) $(wildcard */*.h)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	mkdir -p $(BINDIR)
	$(CXX) $(CCFLAGS) $(WARNFLAGS) $(OPTIMIZERS) $(INCLUDE_DIRS) -o $@ $(OBJS) $(LDFLAGS_DGS) $(LDFLAGS)

crypto_stream_salsa20_amd64_xmm6.o: crypto_stream_salsa20_amd64_xmm6.s
	$(CXX) $(CCFLAGS) -c $?

%.o: %.cpp
	$(CXX) $(CCFLAGS) $(WARNFLAGS) $(OPTIMIZERS) $(INCLUDE_DIRS) -c $< -o $@

clean:
	$(RM) -rf $(TARGET) $(TARGET).dSYM $(OBJS)
