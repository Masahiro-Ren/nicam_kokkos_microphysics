CXX = clang++
INCLUDE_PATH = ../include

CXXFLAGS = --std=c++17 -I$(INCLUDE_PATH)


OBJS = 	\
	problemsize.o \
	data_io.o \
	mod_debug.o \
	mod_satadjust.o \
	main.o \

TARGET = physicskernel_microphysics.exe

$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(TARGET)