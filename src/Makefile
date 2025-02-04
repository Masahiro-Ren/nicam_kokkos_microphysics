KOKKOS_PATH = /home/z42286z/kokkos-4.4.01
# KOKKOS_SRC_PATH = ${KOKKOS_PATH}
# SRC = $(wildcard /home/z42286z/kokkos_stream/*.cpp)
# vpath %.cpp $(sort $(dir $(SRC)))

INCLUDE_PATH = ../include

default: build
	echo "Start Build"

CXX = FCCpx
CXXFLAGS = -Ofast -mcpu=a64fx -fvectorize -ffj-zfill=100 -ffj-prefetch-sequential=soft -ffj-prefetch-line=8 -ffj-prefetch-line-L2=16 -ffj-lst=t -I$(INCLUDE_PATH)
LINK = ${CXX}
LDFLAGS = -Nclang -fopenmp ${CXXFLAGS}
EXE = phy.exe
KOKKOS_DEVICES = "OpenMP"
# KOKKOS_ARCH = "A64FX"

DEPFLAGS = -M

# OBJ = $(notdir $(SRC:.cpp=.o))
OBJ = \
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

LIB =

include $(KOKKOS_PATH)/Makefile.kokkosfx

build: $(EXE)

test: $(EXE)
	./$(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean 
	rm -rf *.o *.cuda $(EXE) *.tmp *.lst desul

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)
