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

#ifndef PDSM_SNF_ALGORITHMS_H
#define PDSM_SNF_ALGORITHMS_H

#include <iostream>
#include <cassert>
#include <boost/optional.hpp>

namespace gyoza {
    namespace algorithms {

        template< typename _Matrix_T >
        // Adopted from Eigen::FullPivLU...
        // _Matrix_T is assumed to be an eigen matrix type.
        class FullPivotSNF {
        public:
            typedef _Matrix_T Matrix_T;
            // Assumed to be an exact field type.
            typedef typename Matrix_T::Scalar Scalar;
            typedef typename Matrix_T::Index Index;

            FullPivotSNF( bool do_inverses = false ) : is_initialized( false ), stores_inverses( do_inverses ) { ; }

            FullPivotSNF( const Matrix_T& matrix, bool do_inverses = false ) : m_snf( matrix.rows(), matrix.cols()),
                                                                               is_initialized( false ),
                                                                               stores_inverses( do_inverses ) {
                original_mat = matrix;
                compute( matrix );
            }

            FullPivotSNF& compute( const Matrix_T& matrix );

            Index rank() const {
                // TODO: upgrade to exception
                if ( !is_initialized ) {
                    throw ("snf solver uninitialized!");
                }
                return nonzero_pivots;
            }

            Index rows() const {
                if ( !is_initialized ) {
                    throw ("snf solver uninitialized!");
                }
                return m_snf.rows();
            }

            Index cols() const {
                if ( !is_initialized ) {
                    throw ("snf solver uninitialized!");
                }
                return m_snf.cols();
            }

            Matrix_T get_mat() const {
                if ( !is_initialized ) {
                    throw ("snf solver uninitialized!");
                }
                return original_mat;
            }

            Matrix_T get_snf() const {
                // TODO: upgrade to exception
                if ( !is_initialized ) {
                    throw ("snf solver uninitialized!");
                }
                return m_snf;
            }

            Matrix_T get_P() const {
                if ( !is_initialized ) {
                    throw ("snf solver uninitialized!");
                }
                return P;
            }

            Matrix_T get_Q() const {
                if ( !is_initialized ) {
                    throw ("snf solver uninitialized!");
                }
                return Q;
            }

            Matrix_T get_P_inv() const {
                if ( !is_initialized || !stores_inverses ) {
                    throw ("snf solver did not compute inverses!");
                }
                return P_inv;
            }

            Matrix_T get_Q_inv() const {
                if ( !is_initialized || !stores_inverses ) {
                    throw ("snf solver did not compute inverses!");
                }
                return Q_inv;
            }

            boost::optional< Matrix_T > get_inverse() const {
                if ( !is_initialized ) {
                    throw ("snf solver uninitialized!");
                }
                if ( rows() != cols()) {
                    throw ("cannot invert nonsquare matrix!");
                }

                if ( nonzero_pivots == rows()) {
                    return (Q * P).eval();
                }
                // not invertible
                return boost::none;
            }

            void print_status( std::ostream& os ) const {
                os << "init: " << is_initialized << "\n";
                os << "inverses: " << stores_inverses << "\n";
            }


        protected:
            Matrix_T original_mat;
            Matrix_T m_snf;
            bool is_initialized;
            bool stores_inverses;

            Index nonzero_pivots;

            Matrix_T P;
            Matrix_T Q;

            Matrix_T P_inv;
            Matrix_T Q_inv;

            // m_snf = P * matrix * Q
        };

    }
}

template< typename Matrix_T >
bool transform_checker( const Matrix_T& matrix,
                        const Matrix_T& matrixprime,
                        const Matrix_T& P,
                        const Matrix_T& Q ) {
    return (P * matrix * Q) == (matrixprime);
}

template< typename Matrix_T >
bool transform_checker_inv( const Matrix_T& matrix,
                            const Matrix_T& matrixprime,
                            const Matrix_T& P_inv,
                            const Matrix_T& Q_inv ) {
    return (matrix) == (P_inv * matrixprime * Q_inv);
}


namespace gyoza {
    namespace algorithms {
        template< typename _Matrix_T >
        FullPivotSNF< _Matrix_T >&
        FullPivotSNF< _Matrix_T >::compute( const Matrix_T& matrix ) {
            original_mat = matrix;
            m_snf = matrix;

            const Index size = matrix.diagonalSize();
            nonzero_pivots = size;
            const Index rows = matrix.rows();
            const Index cols = matrix.cols();

            P = Matrix_T::Identity( rows, rows );
            Q = Matrix_T::Identity( cols, cols );
            if ( stores_inverses ) {
                P_inv = Matrix_T::Identity( rows, rows );
                Q_inv = Matrix_T::Identity( cols, cols );
            }

            for ( Index current_index = 0; current_index < size; ++current_index ) {
                Index row_piv, col_piv;
                Scalar coeff = m_snf.bottomRightCorner( rows - current_index,
                                                        cols - current_index ).maxCoeff( &row_piv, &col_piv );
                row_piv += current_index;
                col_piv += current_index;

                if ( coeff == Scalar( 0 )) {
                    nonzero_pivots = current_index;
                    // zero block left
                    break;
                }

                if ( current_index != row_piv ) {
                    m_snf.row( current_index ).swap( m_snf.row( row_piv ));
                    P.row( current_index ).swap( P.row( row_piv ));
                    if ( stores_inverses ) {
                        P_inv.col( current_index ).swap( P_inv.col( row_piv ));
                    }
                }
                if ( current_index != col_piv ) {
                    m_snf.col( current_index ).swap( m_snf.col( col_piv ));
                    Q.col( current_index ).swap( Q.col( col_piv ));
                    if ( stores_inverses ) {
                        Q_inv.row( current_index ).swap( Q_inv.row( col_piv ));
                    }
                }

                // gaussian elimination in both directions.
                Index rsize = rows - current_index;
                Index csize = cols - current_index;


                Scalar mult = m_snf.coeff( current_index, current_index );
                m_snf.col( current_index ) /= mult;
                Q.col( current_index ) /= mult;
                if ( stores_inverses ) {
                    Q_inv.row( current_index ) *= mult;
                }

                if ( rsize > 1 ) {
                    // Naive version:
                    // for (Index kk = current_index + 1; kk < rsize; ++kk){
                    //   Scalar rmult = m_snf.coeff(kk, current_index);
                    //   m_snf.row(kk) -= rmult * m_snf.row(current_index);
                    //   P.row(kk) -= rmult * P.row(current_index);
                    //   P.col(current_index) += P_inv.col(kk) * rmult;
                    // }

                    auto coltozero = m_snf.col( current_index ).tail( rsize - 1 );
                    P.block( current_index + 1, 0, rsize - 1, rows ) -=
                            coltozero * P.row( current_index );
                    if ( stores_inverses ) {
                        P_inv.col( current_index ) +=
                                P_inv.block( 0, current_index + 1, rows, rsize - 1 ) * coltozero;
                    }

                    // warning to future developers..
                    // operate on m_snf AFTER on P.
                    // otherwise, coltozero will be changed and the operation on P
                    // above will no longer be correct.
                    m_snf.block( current_index + 1, current_index, rsize - 1, csize ) -=
                            coltozero * m_snf.row( current_index ).tail( csize );
                }

                if ( csize > 1 ) {
                    // for (Index kk = current_index + 1; kk < csize; ++kk) {
                    //   Scalar cmult = m_snf.coeff(current_index,kk);
                    //   m_snf.col(kk) -= cmult * m_snf.col(current_index);
                    //   Q.col(kk) -= cmult * Q.col(current_index);
                    //   Q_inv.row(current_index) += Q_inv.row(kk) * cmult;
                    // }

                    auto rowtozero = m_snf.row( current_index ).tail( csize - 1 );
                    Q.block( 0, current_index + 1, cols, csize - 1 ) -=
                            Q.col( current_index ) * rowtozero;

                    if ( stores_inverses ) {
                        Q_inv.row( current_index ) += rowtozero *
                                                      Q_inv.block( current_index + 1, 0, csize - 1, cols );
                    }

                    // Same warning as above.
                    m_snf.row( current_index ).tail( csize - 1 ).setZero();
                }

#ifdef DEBUG
                assert( transform_checker(matrix, m_snf,P,Q) );
                if (stores_inverses) {
                  assert( transform_checker_inv(matrix, m_snf,P_inv,Q_inv) );
                }
#endif
            }
            is_initialized = true;
            return *this;
        }

    }
}

#endif //PDSM_SNF_ALGORITHMS_H
