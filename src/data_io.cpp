#include "data_io.h"

namespace DATA_IO{

char BUF_1D[SIZE_BUF_1D];
char BUF_2D[SIZE_BUF_2D];

char BUF_3D[SIZE_BUF_3D];
char BUF_3D2[SIZE_BUF_3D2];

char BUF_4D[SIZE_BUF_4D];
char BUF_4D2[SIZE_BUF_4D2];

void read_data_1d(const std::string& filename, double arr1d[ADM_kall])
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error opening " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    if(!infile.read(BUF_1D, SIZE_BUF_1D))
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    double* cast_array = reinterpret_cast<double*>(BUF_1D);

    for(int i = 0; i < ADM_kall; i++)
        arr1d[i] = cast_array[i];
}

void read_data_2d(const std::string& filename, double arr2d[ADM_lall][ADM_gall_in])
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error opening " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    if(!infile.read(BUF_2D, SIZE_BUF_2D))
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    double* cast_array = reinterpret_cast<double*>(BUF_2D);

    for(int j = 0; j < ADM_lall; j++)
    {

        for(int i = 0; i < ADM_gall_in; i++)
        {
            arr2d[j][i] = cast_array[j * ADM_gall_in + i];
        }
    }
}

void read_data_3d(const std::string& filename, double arr3d[ADM_lall][ADM_kall][ADM_gall_in])
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error opening " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    if(!infile.read(BUF_3D, SIZE_BUF_3D))
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    double* cast_array = reinterpret_cast<double*>(BUF_3D);

    size_t IJ = ADM_kall * ADM_gall_in;

    for(int k = 0; k < ADM_lall; k++)
    {
        for(int j = 0; j < ADM_kall; j++)
        {
            for(int i = 0; i < ADM_gall_in; i++)
            {
                arr3d[k][j][i] = cast_array[ k * IJ + j * ADM_gall_in + i ];
            }
        }
    }
}

void read_data_3d(const std::string& filename, double arr3d[ADM_lall][ADM_KNONE][ADM_gall_in])
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error opening " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    if(!infile.read(BUF_3D2, SIZE_BUF_3D2))
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    double* cast_array = reinterpret_cast<double*>(BUF_3D2);

    size_t IJ = ADM_KNONE * ADM_gall_in;

    for(int k = 0; k < ADM_lall; k++)
    {
        for(int j = 0; j < ADM_KNONE; j++)
        {
            for(int i = 0; i < ADM_gall_in; i++)
            {
                arr3d[k][j][i] = cast_array[ k * IJ + j * ADM_gall_in + i ];
            }
        }
    }
}

void read_data_4d(const std::string& filename, double arr4d[ADM_lall][TRC_VMAX][ADM_kall][ADM_gall_in])
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error opening " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    if(!infile.read(BUF_4D, SIZE_BUF_4D))
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    double* cast_array = reinterpret_cast<double*>(BUF_4D);

    size_t IJK = TRC_VMAX * ADM_kall * ADM_gall_in;
    size_t IJ = ADM_kall * ADM_gall_in;

    for(int l = 0; l < ADM_lall; l++)
    {
        for(int k = 0; k < TRC_VMAX; k++)
        {
            for(int j = 0; j < ADM_kall; j++)
            {
                for(int i = 0; i < ADM_gall_in; i++)
                {
                    arr4d[l][k][j][i] = cast_array[l * IJK + k * IJ + j * ADM_gall_in + i];
                }
            }
        }
    }
}

void read_data_4d(const std::string& filename, double arr4d[2][ADM_lall][ADM_KNONE][ADM_gall_in])
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error opening " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    if(!infile.read(BUF_4D2, SIZE_BUF_4D2))
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    double* cast_array = reinterpret_cast<double*>(BUF_4D2);

    size_t IJK = ADM_lall * ADM_KNONE * ADM_gall_in;
    size_t IJ = ADM_KNONE * ADM_gall_in;

    for(int l = 0; l < 2; l++)
    {
        for(int k = 0; k < ADM_lall; k++)
        {
            for(int j = 0; j < ADM_KNONE; j++)
            {
                for(int i = 0; i < ADM_gall_in; i++)
                {
                    arr4d[l][k][j][i] = cast_array[l * IJK + k * IJ + j * ADM_gall_in + i];
                }
            }
        }
    }
}

};
