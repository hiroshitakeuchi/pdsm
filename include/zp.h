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

#ifndef PDSM_ZP_H
#define PDSM_ZP_H

#include <iostream>
#include <vector>
#include <cmath>

namespace Core {
    // Zp is the number class of integers modulo p.
    // implemented by int.
    constexpr int MODULO = 1009;

    // return gcd s.t. ax + by = gcd(a, b)
    int euclid( int a, int b, int& x, int& y );

    int inverse_Zp( int num );

    std::vector< int > return_inv_Zp();

    inline const std::vector< int > INVERSE_ZP = return_inv_Zp();

    class Zp {
    public:
        Zp() { num = 0; };

        Zp( int x ) {
            num = ((x % MODULO) + MODULO) % MODULO; // unified but slow?
        };

        Zp( const Zp& ) = default;

        Zp( Zp&& ) = default;

        Zp& operator=( const Zp& ) = default;

        Zp& operator=( Zp&& ) = default;

        Zp& operator=( int other );

        Zp& operator+=( const Zp& other );

        Zp& operator-=( const Zp& other );

        Zp& operator*=( const Zp& other );

        Zp& operator/=( const Zp& other );

        bool operator==( const Zp& other ) const;

        bool operator!=( const Zp& other ) const;

        bool operator>( const Zp& other ) const;

        bool operator<( const Zp& other ) const;

        bool operator>=( const Zp& other ) const;

        bool operator<=( const Zp& other ) const;

        explicit operator int() const;

        friend std::ostream& operator<<( std::ostream& os, const Zp& that );

        friend void operator>>( std::istream& is, Zp& that );

        int num;
    };


    inline Zp operator%( const Zp& lhs, const Zp& rhs ) {
        Zp temp( lhs.num % rhs.num );
        return temp;
    }

    inline Zp operator+( Zp lhs, const Zp& rhs ) {
        lhs += rhs;
        return lhs;
    }

    inline Zp operator-( Zp lhs, const Zp& rhs ) {
        lhs -= rhs;
        return lhs;
    }

    inline Zp operator-( const Zp& rhs ) {
        Zp lhs((rhs.num * (-1) + MODULO) % MODULO );
        return lhs;
    }

    inline Zp operator*( Zp lhs, const Zp& rhs ) {
        lhs *= rhs;
        return lhs;
    }

    inline Zp operator/( Zp lhs, const Zp& rhs ) {
        lhs /= rhs;
        return lhs;
    }

    inline Zp abs( Zp rhs ) {
        return rhs;
    }

    inline int euclid( int a, int b, int& x, int& y ) {
        if ( b == 0 ) {
            x = 1;
            y = 0;
            return a;
        }
        int gcd = euclid( b, a % b, y, x );
        y -= a / b * x;
        return gcd;
    }

    inline int inverse_Zp( int num ) {
        int x, y;
        int gcd = euclid( MODULO, num, x, y );
        return ((y % MODULO) + MODULO) % MODULO;
    }

    inline std::vector< int > return_inv_Zp() {
        std::vector< int > inv_vec( MODULO );
        for ( int i = 0; i < MODULO; ++i ) {
            inv_vec.at( i ) = inverse_Zp( i );
        }
        return inv_vec;
    }

    inline Zp& Zp::operator=( int other ) {
        num = other % MODULO;
        return *this;
    }

    inline Zp& Zp::operator+=( const Zp& other ) {
        num = (num + other.num) % MODULO;
        return *this;
    }

    inline Zp& Zp::operator-=( const Zp& other ) {
        num = (num - other.num + MODULO) % MODULO;
        return *this;
    }

    inline Zp& Zp::operator*=( const Zp& other ) {
        num = (num * other.num) % MODULO;
        return *this;
    }

    inline Zp& Zp::operator/=( const Zp& other ) {
        num = (num * INVERSE_ZP[other.num]) % MODULO;
        return *this;
    }

    inline bool Zp::operator==( const Zp& other ) const {
        return (num == other.num);
    }

    inline bool Zp::operator!=( const Zp& other ) const {
        return !operator==( other );
    }

    inline bool Zp::operator>( const Zp& other ) const {
        return (num > other.num);
    }

    inline bool Zp::operator<( const Zp& other ) const {
        return (num < other.num);
    }

    inline bool Zp::operator>=( const Zp& other ) const {
        return (num >= other.num);
    }

    inline bool Zp::operator<=( const Zp& other ) const {
        return (num <= other.num);
    }

    inline Zp::operator int() const {
        return num % MODULO;
    }


    inline std::ostream& operator<<( std::ostream& os, const Zp& that ) {
        os << that.num;
        return os;
    }

    inline void operator>>( std::istream& is, Zp& that ) {
        is >> that.num;
        that.num = ((that.num % MODULO) + MODULO) % MODULO;
    }

}

#endif //PDSM_ZP_H
