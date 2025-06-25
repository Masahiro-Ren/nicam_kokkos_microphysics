if(ARCHITECTURE STREQUAL "A64FX")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Nclang -Ofast -mcpu=a64fx -fvectorize \ 
#                        -ffj-zfill=100 -ffj-prefetch-sequential=soft \
#                        -ffj-prefetch-line=8 -ffj-prefetch-line-L2=16 \
#                        -ffj-lst=t \
#                        -DUSEOPENMP")
    target_compile_options(phy.exe PRIVATE
                            -Nclang -Ofast -mcpu=a64fx -fvectorize
                            -ffj-zfill=100 -ffj-prefetch-sequential=soft 
                            -ffj-prefetch-line=8 -ffj-prefetch-line-L2=16 
                            -ffj-lst=t -DUSEOPENMP
                            )
elseif(ARCHITECTURE STREQUAL "INTEL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize \ 
                        -march=cascadelake \
                        -DUSEOPENMP")
elseif(ARCHITECTURE STREQUAL "CUDA")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()