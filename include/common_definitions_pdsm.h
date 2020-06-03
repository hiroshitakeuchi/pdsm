/*  Copyright (C) 2018 Hiroshi Takeuchi

    This file is part of pdsm.

    pdsm is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see https://www.gnu.org/licenses/
*/

#ifndef PDSM_COMMON_DEFINITIONS_PDSM_H
#define PDSM_COMMON_DEFINITIONS_PDSM_H

#include <iostream>
#include <vector>
#include <fstream>
#include <functional>
#include <cmath>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <cassert>
#include <set>
#include <utility>
#include <float.h>
#include <boost/progress.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include "zp.h"
#include "zptraits.h"
#include "common_definitions_zp.h"
#include "snf_algorithms.h"

namespace Pdsm {
    typedef int Index;
    typedef int16_t Dimension;
    typedef std::vector< Index > IndexColumn;
    typedef double Radius;
    typedef std::vector< Radius > Radii;
    typedef std::vector< std::pair< Pdsm::Index, Pdsm::Radius >> CriticalIndexRadiusPairs;
    typedef gyoza::ZpSparseMatrix SparseMatrix;
    typedef gyoza::ZpSparseVector SparseColumn;
    typedef gyoza::ZpMatrix DenseMatrix;
//    typedef gyoza::ZpMatrix DenseColumn;

    // max Index whose Radius is lower than value
    Index max_index_under_value( const Radii& rad, const double value ) {
        size_t size = rad.size();
        int i;
        for ( i = 0; i < size; ++i ) {
            if ( rad.at( i ) > value ) {
                break;
            }
        }
        return i - 1;
    }

    void sort_and_unique( Radii& rad ) {
        std::sort( rad.begin(), rad.end());
        rad.erase( std::unique( rad.begin(), rad.end()), rad.end());
    }

    template< typename STLofIndex >
    void print_STL( STLofIndex set ) {
        Index idx = 0;
        for ( const auto& a : set ) {
            std::cout << a << ", ";
            ++idx;
        }
        if ( idx != 0 ) {
            std::cout << "\b\b" << std::endl;
        } else {
            std::cout << std::endl;
        }
    }

}

#endif //PDSM_COMMON_DEFINITIONS_PDSM_H
