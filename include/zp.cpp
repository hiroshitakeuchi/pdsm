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

#include "zp.h"
#include <cmath>

namespace Core {
    int euclid( int a, int b, int& x, int& y ) {
        // a<=b
        if (a == 0) {
            x = 0;
            y = 1;
            return b;
        }
        int gcd = euclid(b%a, a, y, x);
        // y,x in former step (b%a*y + a*(x+b/a*y) = gcd)
        x -= b/a * y;
        return gcd;
    }

    int inverse_Zp(int num){
        int x, y;
        int gcd = euclid(MODULO, num, x, y); // gcd = 1
        // 1 = MODULO*x + num*y = num*y (mod MODULO)
        return ((y % MODULO) + MODULO) % MODULO;
    }

    std::vector<int> return_inv_Zp(){
        std::vector<int> inv_vec(MODULO);
        for ( int i = 0; i < MODULO ; ++i ) {
            inv_vec.at(i) = inverse_Zp(i);
        }
        return inv_vec;
    }

    Zp& Zp::operator=( int other ) {
        num = other % MODULO;
        return *this;
    }

    Zp& Zp::operator+=( const Zp& other ) {
        num = (num + other.num) % MODULO;
        return *this;
    }

    Zp& Zp::operator-=( const Zp& other ) {
        num = (num - other.num + MODULO) % MODULO;
        return *this;
    }

    Zp& Zp::operator*=( const Zp& other ) {
        num = (num * other.num) % MODULO;
        return *this;
    }

    Zp& Zp::operator/=( const Zp& other ) {
        num = (num * INVERSE_ZP[other.num]) % MODULO;
        return *this;
    }

    bool Zp::operator==( const Zp& other ) const {
        return (num == other.num);
    }

    bool Zp::operator!=( const Zp& other ) const {
        return !operator==( other );
    }

    bool Zp::operator>( const Zp& other ) const {
        return (num > other.num);
    }

    bool Zp::operator<( const Zp& other ) const {
        return (num < other.num);
    }

    bool Zp::operator>=( const Zp& other ) const {
        return (num >= other.num);
    }

    bool Zp::operator<=( const Zp& other ) const {
        return (num <= other.num);
    }

    Zp::operator int() const {
        return num % MODULO;
    }


    std::ostream& operator<<( std::ostream& os, const Zp& that ) {
        os << that.num;
        return os;
    }

    void operator>>( std::istream& is, Zp& that ) {
        is >> that.num;
        that.num = ((that.num % MODULO) + MODULO) % MODULO;
    }
}