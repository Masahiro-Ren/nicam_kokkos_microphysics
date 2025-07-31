if(ARCHITECTURE STREQUAL "A64FX")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Nclang -Ofast -mcpu=a64fx -fvectorize \ 
#                        -ffj-zfill=100 -ffj-prefetch-sequential=soft \
#                        -ffj-prefetch-line=8 -ffj-prefetch-line-L2=16 \
#                        -ffj-lst=t \
#                        -DUSEOPENMP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Nclang -Ofast -mcpu=a64fx -fvectorize -ffj-lst=t -DUSE_OPENMP")
elseif(ARCHITECTURE STREQUAL "INTEL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -march=cascadelake -DUSE_OPENMP")
elseif(ARCHITECTURE STREQUAL "CUDA")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DUSE_CUDA")
elseif(ARCHITECTURE STREQUAL "GRACE_HOPPER")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -mp -DUSE_CUDA")
endif()