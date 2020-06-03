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

#ifndef PDSM_ZPTRAITS_H
#define PDSM_ZPTRAITS_H

#include "zp.h"
#include <Eigen/Core>

namespace Eigen {
    template<> struct NumTraits<Core::Zp> {
        typedef Core::Zp Real;
        typedef Core::Zp Nested;
        typedef Core::Zp Literal;
        enum {
            IsComplex = 0,
            IsInteger = 1,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 3,
            MulCost = 3,
        };
        static inline Real dummy_precision() { return Real(0); }
        static inline int digits10() { return 1; } // return internal::default_digits10_impl<T>::run(); }

    };

}

#endif //PDSM_ZPTRAITS_H
