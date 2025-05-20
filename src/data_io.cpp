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
    {
        arr1d[i] = cast_array[i];
    }
}

/**
 * Kokkos ver.
 */
void read_data_1d(const std::string& filename, View<double*>& arr1d)
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

    size_t dim0 = arr1d.extent_int(0);

    // RangePolicy<Schedule<Kokkos::Static>, EXE_SPACE> range_1d(0, dim0);
    Kokkos::parallel_for("read_data_1d", RangePolicy<>(0,dim0), 
    KOKKOS_LAMBDA(const size_t i){
        arr1d(i) = cast_array[i];
    });
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

/**
 * Kokkos ver.
 */
void read_data_2d(const std::string& filename, View<double**>& arr2d)
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

    size_t dim0 = arr2d.extent_int(0);
    size_t dim1 = arr2d.extent_int(1);

    // MDRangePolicy<EXE_SPACE, Kokkos::Rank<2>> range_2d({0,0},{dim0, dim1});
    Kokkos::parallel_for("read_data_2d", MDRangePolicy<Kokkos::Rank<2>>({0,0},{dim0,dim1}), 
    KOKKOS_LAMBDA(const size_t j, const size_t i){
        arr2d(j, i) = cast_array[j * ADM_gall_in + i];
    });
    
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

/**
 * Kokkos ver.
 */
void read_data_3d(const std::string& filename, View<double***>& arr3d)
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error opening " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    size_t dim0 = arr3d.extent_int(0);
    size_t dim1 = arr3d.extent_int(1);
    size_t dim2 = arr3d.extent_int(2);

    size_t IJ = dim1 * dim2;

    double* cast_array = nullptr;

    /* May have better implementation*/
    if(dim1 == ADM_kall)
    {
        if(!infile.read(BUF_3D, SIZE_BUF_3D))
        {
            std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
            DEBUG::ADM_Proc_stop();
        }

        cast_array = reinterpret_cast<double*>(BUF_3D);
    }
    else if(dim1 == ADM_KNONE)
    {
        if(!infile.read(BUF_3D2, SIZE_BUF_3D2))
        {
            std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
            DEBUG::ADM_Proc_stop();
        }

        cast_array = reinterpret_cast<double*>(BUF_3D2);
    }

    if(cast_array == nullptr)
    {
        std::cerr << __PRETTY_FUNCTION__ << " cast_array is nullptr while reading " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    // MDRangePolicy<EXE_SPACE, Kokkos::Rank<3>> range_3d({0,0,0}, {dim0, dim1, dim2}); 
    Kokkos::parallel_for("read_data_3d", MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{dim0,dim1,dim2}), 
    KOKKOS_LAMBDA(const size_t k, const size_t j, const size_t i){
        arr3d(k, j, i) = cast_array[k * IJ + j * ADM_gall_in + i];
    });
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

/**
 * Kokkos ver.
 */
void read_data_4d(const std::string& filename, View<double****>& arr4d)
{
    std::ifstream infile(filename, std::ios::binary);

    if(!infile)
    {
        std::cerr << __PRETTY_FUNCTION__ << " Error opening " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    size_t dim0 = arr4d.extent_int(0);
    size_t dim1 = arr4d.extent_int(1);
    size_t dim2 = arr4d.extent_int(2);
    size_t dim3 = arr4d.extent_int(3);

    size_t IJK = dim1 * dim2 * dim3;
    size_t IJ = dim2 * dim3;

    double* cast_array = nullptr;

    if(dim2 == ADM_kall)
    {
        if(!infile.read(BUF_4D, SIZE_BUF_4D))
        {
            std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
            DEBUG::ADM_Proc_stop();
        }

        cast_array = reinterpret_cast<double*>(BUF_4D);
    }
    else if(dim2 == ADM_KNONE)
    {
        if(!infile.read(BUF_4D2, SIZE_BUF_4D2))
        {
            std::cerr << __PRETTY_FUNCTION__ << " Error reading " << filename << std::endl;
            DEBUG::ADM_Proc_stop();
        }

        cast_array = reinterpret_cast<double*>(BUF_4D2);
    }

    if(cast_array == nullptr)
    {
        std::cerr << __PRETTY_FUNCTION__ << " cast_array is nullptr while reading " << filename << std::endl;
        DEBUG::ADM_Proc_stop();
    }

    // MDRangePolicy<EXE_SPACE, Kokkos::Rank<4>> range_4d({0,0,0,0}, {dim0, dim1, dim2, dim3}); 
    Kokkos::parallel_for("read_data_4d", MDRangePolicy<Kokkos::Rank<4>>({0,0,0,0},{dim0,dim1,dim2,dim3}), 
    KOKKOS_LAMBDA(const size_t l, const size_t k, const size_t j, const size_t i){
        arr4d(l, k, j, i) = cast_array[l * IJK + k * IJ + j * ADM_gall_in + i];
    });
}

}
