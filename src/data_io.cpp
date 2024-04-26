#include "data_io.h"

void Data_IO::read(VAR_LIST var_name, Vec2d<double>& data)
{
    std::string prefix = "./data/";
    std::string suffix = ".dat";
    std::string file_name = prefix + var_names[var_name] + suffix;

    std::cout << "Reading data from " << file_name << std::endl;
}

void Data_IO::close(VAR_LIST var_name)
{
    file_list[var_name].close();
}

void Data_IO::closeall()
{
    for(auto& file : file_list)
    {
        file.close();
    }
}