# CXX = g++
CXX = FCCpx
INCLUDE_PATH = ../include

# FLAGS for a64fx
A64FX_FLAGS = --std=c++17 -g -O0

DEBUG_FLAGS = --std=c++17 -g -O0 -Wall -Wextra -Wpedantic -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=float-divide-by-zero,float-cast-overflow

CXXFLAGS = $(A64FX_FLAGS) -I$(INCLUDE_PATH) 
# CXXFLAGS = $(DEBUG_FLAGS) -I$(INCLUDE_PATH)

OBJS = 	\
	problemsize.o          \
	data_io.o              \
	mod_debug.o            \
	mod_vadv1d.o           \
	mod_thrmdyn.o          \
	mod_satadjust.o        \
	mod_precip_transport.o \
	mod_mp_nsw6.o          \
	mod_mp_driver.o        \
	main.o

TARGET = physicskernel_microphysics.exe

$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(TARGET)