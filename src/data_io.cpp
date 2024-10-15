#include "data_io.h"

namespace DATA_IO{

char BUF_1D[SIZE_BUF_1D];
char BUF_2D[SIZE_BUF_2D];

char BUF_3D[SIZE_BUF_3D];
char BUF_3D[SIZE_BUF_3D2];

char BUF_4D[SIZE_BUF_4D];
char BUF_4D[SIZE_BUF_4D2];

void read_data_1d(const std::string& filename, double arr1d[ADM_kall])
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << "Error opeing file. \n";
        DEBUG::ADM_Proc_stop();
    }

    if(!infile.read(BUF_1D, SIZE_BUF_1D))
    {
        std::cerr << "Error reading file. \n";
        DEBUG::ADM_Proc_stop();
    }

    double* cast_array = reinterpret_cast<double*>(BUF_1D);

    for(int i = 0; i < ADM_kall; i++)
        arr1d[i] = cast_array[i];
}

void read_data_2d(const std::string& filename, double arr2d[ADM_kall][ADM_gall_in])
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << "Error opeing file. \n";
        DEBUG::ADM_Proc_stop();
    }

    if(!infile.read(BUF_2D, SIZE_BUF_2D))
    {
        std::cerr << "Error reading file. \n";
        DEBUG::ADM_Proc_stop();
    }

    double* cast_array = reinterpret_cast<double*>(BUF_2D);

    for(int j = 0; j < ADM_kall; j++)
    {

        for(int i = 0; i < ADM_gall_in; i++)
            arr2d[j][i] = cast_array[j * ADM_gall_in + i];
    }
}
};
