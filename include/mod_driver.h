#pragma once

#include "problemsize.h"
#include "mod_debug.h"

using namespace PROBLEM_SIZE;
using namespace DEBUG;

namespace DRIVER{

	void mp_init(const std::string& MP_TYPE_in);

	template<size_t ijdim, size_t kdim>
	void mp_dirver();

}