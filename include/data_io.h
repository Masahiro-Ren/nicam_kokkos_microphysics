#ifndef DATA_IO
#define DATA_IO

#include "problemsize.h"
#include <fstream>

enum VAR_LIST {
    RHOG,
    RHOGVX,
    RHOGVY,
    RHOGVZ,
    RHOGW,
    RHOGE,
    RHOGQ_LSWP,
    VX,
    VY,
    VZ,
    W,
    UNCCN,
    RHO,
    TEM,
    PRE,
    Q_LSWP,
    PRECIP_MP,
    PRECIP1_MP,
    PRECIP2_MP,
    RHOEIN_PRECIP_MP,
    LH_PRECIP_MP,
    RHOPHI_PRECIP_MP,
    RHOKIN_PRECIP_MP,
    FRHOGE_AF,
    FRHOGQV_AF,
    FRHOGE_RAD,
    QKE,
    GSGAM2,
    GSGAM2H,
    GAM2,
    GAM2H,
    IX,
    IY,
    IZ,
    JX,
    JY,
    JZ,
    Z,
    ZH,
    GPREC,
    CBMFX
};

class Data_IO {
private:
    std::vector<std::string> var_names;
    std::vector<std::ifstream> file_list;
public:
    Data_IO()
    {
        var_names = {
                    "rhog",
                    "rhogvx",
                    "rhogvy",
                    "rhogvz",
                    "rhogw",
                    "rhoge",
                    "rhogq_Lswp",
                    "vx",
                    "vy",
                    "vz",
                    "w",
                    "unccn",
                    "rho",
                    "pre",
                    "tem",
                    "q_Lswp",
                    "precip_mp",
                    "precip1_mp",
                    "precip2_mp",
                    "rhoein_precip_mp",
                    "lh_precip_mp",
                    "rhophi_precip_mp",
                    "rhokin_precip_mp",
                    "frhoge_af",
                    "frhoqv_af",
                    "frhoge_rad",
                    "qke",
                    "gsgam2",
                    "gsgam2h",
                    "gam2",
                    "gam2h",
                    "ix",
                    "iy",
                    "iz",
                    "jx",
                    "jy",
                    "jz",
                    "z",
                    "zh",
                    "GPREC",
                    "CBMFX",
                    "qd",
                    "rceff",
                    "rceff_solid",
                    "rceff_cld",
                    "rctop",
                    "rwtop",
                    "tctop",
                    "GDCLW",
                    "GDCFRC"};
    }

    ~Data_IO()
    {
        closeall();
    }

    void read(VAR_LIST var_name, Vec2d<double>& data);
    void close(VAR_LIST var_name);
    void closeall();
};

#endif

