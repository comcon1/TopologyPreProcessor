/** \file global.hpp
 *
 *	\brief This file is used to interact with AutoMake toolchain. It should be
 *	incluided whenever conditional compilation is used.
 *
 *	Well, a relict of  its former glory, isnt't?
 */

#ifndef TPP_GLOBAL_H
#define TPP_GLOBAL_H

#ifdef HAVE_CONFIG_H
#include "../config.h"
#else
#error !! YOU SHOULD RUN CONFIGURE SCRIPT !!
#endif

#include <ostream>
#include <vector>


template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];
        if (i != last) 
            out << ", ";
    }
    out << "]";
    return out;
}


#endif
