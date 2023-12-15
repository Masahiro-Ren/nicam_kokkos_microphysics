#ifndef NDARRAY_H
#define NDARRAY_H

#include <vector>

template <typename T>
using Vec1D = std::vector<T>;

template <typename T>
using Vec2D = std::vector<std::vector<T>>;

template <typename T>
using Vec3D = std::vector<std::vector<std::vector<T>>>;

template <typename T>
using Vec4D = std::vector<std::vector<std::vector<std::vector<T>>>>;

#endif