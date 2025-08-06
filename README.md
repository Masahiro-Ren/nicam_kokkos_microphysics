# Kokkos-Based NICAM Physics Kernel: Microphysics Module

This repository contains the Kokkos-accelerated version of the NICAM Physics Kernel for the microphysics module.

**Porting Path:**  
Fortran → Pure C++ → Kokkos-based Parallel Implementation

**Authorization:**  
Parts of this code are derived from the original [NICAM dckernel 2016](https://github.com/hisashiyashiro/nicam_dckernel_2016), licensed under the BSD 2-Clause license. Permission to redistribute and modify this ported version was granted by the original author.

---

## Requirements

- **CMake** ≥ 3.20  
- **C++17-Compatible Compilers**:  
  - `clang++`, `g++`, `nvcc`, `nvc++`
  - `FCCpx` in clang mode
- **Kokkos** ≥ 4.6.1  
  Installation guide available here: [Kokkos Documentation](https://kokkos.org/kokkos-core-wiki/get-started.html)

---

## Tested Environments

- ✅ **Flow Type I Subsystem** — *Fujitsu A64FX @ Nagoya University*  
- ✅ **Flow Type II Subsystem** — *Intel Cascade Lake + NVIDIA V100 @ Nagoya University*  
- ✅ **Miyabi-G** — *NVIDIA GH200 Superchip @ University of Tokyo*

---

## Build Instructions

### ▶️ Flow Type I (Fujitsu A64FX)

```bash
module load cmake
cmake -B build \
  -DCMAKE_CXX_COMPILER=FCCpx \
  -DKokkos_ROOT=/path/to/your/kokkos/lib64/cmake/Kokkos/ \
  -DARCHITECTURE=A64FX
```

### ▶️ Flow Type II (Intel Cascade Lake + NVIDIA V100)
``` bash
module load gcc/11.3.0 cuda/12.4.1 cmake
cmake -B build \
  -DKokkos_ROOT=/path/to/your/kokkos/lib64/cmake/Kokkos \
  -DARCHITECTURE=V100
```

### ▶️ Miyabi-G (NVIDIA GH200 Superchip)
``` bash
cmake -B build \
  -DCMAKE_CXX_COMPILER=nvc++ \
  -DKokkos_ROOT=/path/to/your/kokkos/lib64/cmake/Kokkos \
  -DARCHITECTURE=GRACE_HOPPER
```

### 🔧 Notes on OpenMP-Only Builds
If you plan to run the program using OpenMP-only, it is recommended to maintain separate Kokkos installations:
``` bash
/kokkos
├── openmp_version      # Built with only OpenMP backend
└── full_version        # Built with all available backends
```
Make sure to point Kokkos_ROOT to the appropriate version during configuration.

## Running Instructions

Copy the executable file to the `run` directory, and create a job script suitable for your own system.

The structure of the `run` directory should be as follows:

```bash
run
├── data             # contains input data
│   └── vgrid        # contains vgrid data
└── ref_verify       # contains verification data
```
**NOTE**: A proper location for storing all input data has not yet been determined. Currently, only the vgrid data and verification data are available in the repository.

## Results Verification
🔹 Results from Original Fortran Code
```
 ### Check ###
 +check  [check_rhog      ] max=  5.0626169922907140E-12,min= -1.4747092436095950E-12,sum= -1.2059661602670250E-10
 +check  [check_rhogvx    ] max=  4.8540726993451240E-11,min= -1.5639045614079800E-11,sum= -5.3378612595696370E-10
 +check  [check_rhogvy    ] max=  2.9647395649590180E-11,min= -3.8473224606150320E-11,sum=  8.7347003859417770E-10
 +check  [check_rhogvz    ] max=  1.7783108319235910E-11,min= -3.0207836232420960E-11,sum=  4.9090170867367460E-10
 +check  [check_rhogw     ] max=  3.6093350530563840E-13,min= -7.6272321791748250E-13,sum=  5.4313126964076830E-12
 +check  [check_rhoge     ] max=  1.4551915228366850E-10,min= -7.1497319893765960E+01,sum= -1.0776720015301340E+07
 +check  [QV              ] max=  4.0932895639233390E-07,min= -2.7807046225052780E-05,sum= -3.7963720189637420E-02
 +check  [QC              ] max=  2.7807046225052680E-05,min= -4.0932905538366300E-07,sum=  3.7965023327871660E-02
 +check  [QR              ] max=  2.0210040839233920E-11,min= -2.1175823681357510E-22,sum=  3.3737971724727460E-09
 +check  [QI              ] max=  1.9067742179908570E-10,min= -5.3637422603355350E-10,sum= -1.0367761562136360E-06
 +check  [QS              ] max=  4.7746961589686530E-11,min= -1.2116527009913440E-10,sum= -2.6910244100958930E-07
 +check  [QG              ] max=  8.2574226106153230E-13,min= -5.5229191280467250E-11,sum= -7.5400921928473070E-10
```
🔹 Results from Kokkos OpenMP Backend (A64FX)
```
Checking Reuslts
Checking [rhog]   Max = 5.0627280145931763e-12; Min = -1.4747092436095954e-12; Sum = -1.2059694909360985e-10;
Checking [rhogvx] Max = 4.8540726993451244e-11; Min = -1.5642598327758606e-11; Sum = -5.3441964425988498e-10;
Checking [rhogvy] Max = 2.9643842935911380e-11; Min = -3.8474112784570025e-11; Sum =  8.7370192989443408e-10;
Checking [rhogvz] Max = 1.7783108319235907e-11; Min = -3.0208724410840659e-11; Sum =  4.9086714951258621e-10;
Checking [rhogw]  Max = 3.6093350530563839e-13; Min = -7.6272321791748254e-13; Sum =  5.4312784440893621e-12;
Checking [rhoge]  Max = 1.4551915228366852e-10; Min = -7.1497319893824169e+01; Sum = -1.0776720015301270e+07;
Checking [QV]     Max = 4.0932895639233385e-07; Min = -2.7807046225070131e-05; Sum = -3.7963720189639810e-02;
Checking [QC]     Max = 2.7807046225862683e-05; Min = -4.0932905538366303e-07; Sum =  3.7965023327889057e-02;
Checking [QR]     Max = 2.0210040839273629e-11; Min = -1.0649436078802570e-15; Sum =  3.3737818754772465e-09;
Checking [QI]     Max = 1.9067742180586193e-10; Min = -5.3637422603355347e-10; Sum = -1.0367761549578148e-06;
Checking [QS]     Max = 4.7746961589739465e-11; Min = -1.2116527009998145e-10; Sum = -2.6910244100750483e-07;
Checking [QG]     Max = 8.2574226106442741e-13; Min = -5.5229191280599598e-11; Sum = -7.5400921928499098e-10;
```
🔹 Results from Kokkos CUDA Backend (V100)
```
Checking Reuslts
Checking [rhog]   Max = 5.0626169922907138e-12; Min = -1.4747092436095954e-12; Sum = -1.2059661602670246e-10;
Checking [rhogvx] Max = 4.8540726993451244e-11; Min = -1.5642598327758606e-11; Sum = -5.3443867070696949e-10;
Checking [rhogvy] Max = 2.9643842935911380e-11; Min = -3.8473224606150325e-11; Sum =  8.7373166305481232e-10;
Checking [rhogvz] Max = 1.7783108319235907e-11; Min = -3.0207836232420959e-11; Sum =  4.9090170867367462e-10;
Checking [rhogw]  Max = 3.6093350530563839e-13; Min = -7.6272321791748254e-13; Sum =  5.4313126964076832e-12;
Checking [rhoge]  Max = 1.1641532182693481e-10; Min = -7.1497319893736858e+01; Sum = -1.0776720015301347e+07;
Checking [QV]     Max = 4.0932895639233385e-07; Min = -2.7807046225050182e-05; Sum = -3.7963720189640206e-02;
Checking [QC]     Max = 2.7807046225049748e-05; Min = -4.0932905538366303e-07; Sum =  3.7965023327872743e-02;
Checking [QR]     Max = 2.0210040839240542e-11; Min = -3.1763735522036263e-22; Sum =  3.3737971724918941e-09;
Checking [QI]     Max = 1.9067742179908567e-10; Min = -5.3637422603185941e-10; Sum = -1.0367761548868546e-06;
Checking [QS]     Max = 4.7746961589633586e-11; Min = -1.2116527009913442e-10; Sum = -2.6910244100756475e-07;
Checking [QG]     Max = 8.2574226107083806e-13; Min = -5.5229191280467249e-11; Sum = -7.5400921928467841e-10;
```