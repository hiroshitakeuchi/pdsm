/*  This file is part of pdsm.
    This file incorporates a modified version of gyoza.
    Modified Date: 2018/10/26
    By Hiroshi Takeuchi

    Copyright (C) 2018 Emerson G. Escolar

    gyoza is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see https://www.gnu.org/licenses/.
*/

#ifndef PDSM_COMMON_DEFINITIONS_ZP_H
#define PDSM_COMMON_DEFINITIONS_ZP_H

#include "zp.h"
#include "zptraits.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace gyoza {
    typedef Eigen::Matrix< Core::Zp, Eigen::Dynamic, Eigen::Dynamic > ZpMatrix;
    typedef Eigen::SparseMatrix< Core::Zp > ZpSparseMatrix;
    typedef Eigen::SparseVector< Core::Zp > ZpSparseVector;
    typedef ZpMatrix::Index Index;
    typedef ZpMatrix::Scalar Scalar;
}

#endif //PDSM_COMMON_DEFINITIONS_ZP_H
