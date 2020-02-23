PROG  := parallelstablesort
PROG2 := makestablesortdata

SRCS  := parallelstablesort.cpp
SRCS2 := makestablesortdata.cpp

OBJS  = parallelstablesort.o
OBJS2 = makestablesortdata.o

DEPS  = parallelstablesort.d
DEPS2 = makestablesortdata.d

VPATH  = src/parallelstablesort src/makestablesortdata
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -mtune=native -march=native -pipe -std=c++17 -fopenmp
LDFLAGS = -L/home/dc1394/oss/tbb/lib/intel64/gcc4.8 -ltbb \
		  -L/home/dc1394/oss/boost_1_72_0/stage/gcc/lib \
		  -lboost_serialization -lboost_system -lboost_thread

all: $(PROG) $(PROG2) ;
#rm -f $(OBJS) $(DEPS)

-include $(DEPS)

$(PROG): $(OBJS)
		$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) -o $@

$(PROG2): $(OBJS2)
		$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -MMD -MP -DDEBUG $<

clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
		rm -f $(PROG2) $(OBJS2) $(DEPS2)
