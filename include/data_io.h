#pragma once

#include "problemsize.h"
#include "mod_debug.h"
#include <fstream>

using namespace PROBLEM_SIZE;

namespace DATA_IO{
// enum VAR_LIST {
//     RHOG,
//     RHOGVX,
//     RHOGVY,
//     RHOGVZ,
//     RHOGW,
//     RHOGE,
//     RHOGQ_LSWP,
//     VX,
//     VY,
//     VZ,
//     W,
//     UNCCN,
//     RHO,
//     TEM,
//     PRE,
//     Q_LSWP,
//     PRECIP_MP,
//     PRECIP1_MP,
//     PRECIP2_MP,
//     RHOEIN_PRECIP_MP,
//     LH_PRECIP_MP,
//     RHOPHI_PRECIP_MP,
//     RHOKIN_PRECIP_MP,
//     FRHOGE_AF,
//     FRHOGQV_AF,
//     FRHOGE_RAD,
//     QKE,
//     GSGAM2,
//     GSGAM2H,
//     GAM2,
//     GAM2H,
//     IX,
//     IY,
//     IZ,
//     JX,
//     JY,
//     JZ,
//     Z,
//     ZH,
//     GPREC,
//     CBMFX
// };
// std::vector<std::string> var_names = {
//                                         "rhog",
//                                         "rhogvx",
//                                         "rhogvy",
//                                         "rhogvz",
//                                         "rhogw",
//                                         "rhoge",
//                                         "rhogq_Lswp",
//                                         "vx",
//                                         "vy",
//                                         "vz",
//                                         "w",
//                                         "unccn",
//                                         "rho",
//                                         "pre",
//                                         "tem",
//                                         "q_Lswp",
//                                         "precip_mp",
//                                         "precip1_mp",
//                                         "precip2_mp",
//                                         "rhoein_precip_mp",
//                                         "lh_precip_mp",
//                                         "rhophi_precip_mp",
//                                         "rhokin_precip_mp",
//                                         "frhoge_af",
//                                         "frhoqv_af",
//                                         "frhoge_rad",
//                                         "qke",
//                                         "gsgam2",
//                                         "gsgam2h",
//                                         "gam2",
//                                         "gam2h",
//                                         "ix",
//                                         "iy",
//                                         "iz",
//                                         "jx",
//                                         "jy",
//                                         "jz",
//                                         "z",
//                                         "zh",
//                                         "GPREC",
//                                         "CBMFX",
//                                         "qd",
//                                         "rceff",
//                                         "rceff_solid",
//                                         "rceff_cld",
//                                         "rctop",
//                                         "rwtop",
//                                         "tctop",
//                                         "GDCLW",
//                                         "GDCFRC"};

constexpr size_t SIZE_BUF_1D = ADM_kall * sizeof(double);
constexpr size_t SIZE_BUF_2D = ADM_lall * ADM_gall_in * sizeof(double);

constexpr size_t SIZE_BUF_3D = ADM_lall * ADM_kall * ADM_gall_in * sizeof(double);
constexpr size_t SIZE_BUF_3D2 = ADM_lall * ADM_KNONE * ADM_gall_in * sizeof(double);

constexpr size_t SIZE_BUF_4D = ADM_lall * TRC_VMAX * ADM_kall * ADM_gall_in * sizeof(double);
constexpr size_t SIZE_BUF_4D2 = 2 * ADM_lall * ADM_KNONE * ADM_gall_in * sizeof(double);

void read_data_1d(const std::string& filename, double arr1d[ADM_kall]);

void read_data_2d(const std::string& filename, double arr2d[ADM_lall][ADM_gall_in]);

void read_data_3d(const std::string& filename, double arr3d[ADM_lall][ADM_kall][ADM_gall_in]);
void read_data_3d(const std::string& filename, double arr3d[ADM_lall][ADM_KNONE][ADM_gall_in]);

void read_data_4d(const std::string& filename, double arr4d[ADM_lall][TRC_VMAX][ADM_kall][ADM_gall_in]);
void read_data_4d(const std::string& filename, double arr4d[2][ADM_lall][ADM_KNONE][ADM_gall_in]);

/**
 * Kokkos ver.
 */
// void read_data_1d(const std::string& filename, View<double*>& arr1d);
void read_data_1d(const std::string& filename, View1D<double, DEFAULT_MEM> arr1d);
void read_data_2d(const std::string& filename, View2D<double, DEFAULT_MEM> arr2d);
void read_data_3d(const std::string& filename, View3D<double, DEFAULT_MEM> arr3d);
void read_data_4d(const std::string& filename, View4D<double, DEFAULT_MEM> arr4d);

}
