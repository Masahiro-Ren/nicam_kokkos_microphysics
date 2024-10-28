# CXX = clang++
CXX = FCCpx
INCLUDE_PATH = ../include

# FLAGS for a64fx
A64FX_FLAGS = -Kfast

CXXFLAGS = --std=c++17 -I$(INCLUDE_PATH) $(A64FX_FLAGS)
# CXXFLAGS = --std=c++17 -I$(INCLUDE_PATH)

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