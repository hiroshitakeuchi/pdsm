/*  This file is part of pdsm.
    This file incorporates a modified version of PHAT.
    Modified Date: 2018/10/26
    By Hiroshi Takeuchi

    Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Michael Kerber, Jan Reininghaus

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PDSM_PERSISTENT_HOMOLOGY_H
#define PDSM_PERSISTENT_HOMOLOGY_H

#include "common_definitions_pdsm.h"

namespace Pdsm{
    // to preserve set_vertex_map
    class IndexMatrix {

    private:
        std::vector< IndexColumn > matrix;

    public:
        Index cols() const {
            return (Index) matrix.size();
        }

        void resize_cols( Index cols ) {
            matrix.resize( cols );
        }

        void get_col( Index idx, IndexColumn& col ) const {
            col = matrix[idx];
        }

        void set_col( Index idx, const IndexColumn& col ) {
            matrix[idx] = col;
        }

        bool is_empty( Index idx ) const {
            return matrix[idx].empty();
        }

        // print the matrix
        void print_matrix() {
            const Index size = this->cols();
            for ( Index idx = 0; idx != size; ++idx ) {
                print_STL< IndexColumn >( matrix[idx] );
            }
        }

        // Remark: call immediately after generating BoundaryMatrix
        void get_simplicial_map_vertex( const IndexColumn& map_vector, IndexMatrix& simplicial_map_vertex ) const {
            const Index size = this->cols();
            simplicial_map_vertex.resize_cols( size );
            for ( Index idx = 0; idx != size; ++idx ) {
                IndexColumn temp_col;
                temp_col.reserve( matrix.size());
                for ( const auto& vertex : matrix[idx] ) {
                    temp_col.push_back( map_vector.at( vertex ));
                }

                // DO NOT SORT (for sgn function)
                simplicial_map_vertex.set_col( idx, temp_col );
            }
        }
    };


    // boundary matrix with coefficient Z/pZ
    class BoundaryMatrix {

    protected:
        std::vector< Dimension > dims;
        SparseMatrix matrix;

    public:
        Index cols() const {
            return (Index) matrix.cols();
        }

        void resize_cols( Index cols ) {
            dims.resize( cols );
            matrix.resize( cols, cols );
        }

        Dimension get_dim( Index idx ) const {
            return dims[idx];
        }

        void set_dim( Index idx, Dimension dim ) {
            dims[idx] = dim;
        }

        void get_col( Index idx, SparseColumn& col ) const {
            col = matrix.col( idx );
        }

        void set_col( Index idx, const SparseColumn& col ) {
            matrix.col( idx ) = col;
            matrix.prune( 0, 0 );
        }

        bool is_zero( Index idx ) {
            matrix.prune( 0, 0 );
            return (matrix.col( idx ).nonZeros() == 0);
        }

        // max nonzero row Index at idx-column
        Index get_max_index( Index idx ) {
            matrix.prune( 0, 0 );
            if ( matrix.col( idx ).nonZeros() == 0 ) {
                return -1;
            } else {
                gyoza::ZpSparseMatrix::ReverseInnerIterator it( matrix, idx );
                return it.row();
            }
        }

        // reduce the max row of SparseColumn 'target' by using SparseColumn 'source'
        // This is used only for the reduction of boundary matrices.
        // Assume1: the max row of source and target should be same.
        // Assume2: the coeff of the max row of source should be 1.
        void reduce( Index source, Index target ) {
            SparseMatrix::ReverseInnerIterator it( matrix, target );
            if ( it ) {
                matrix.col( target ) -= matrix.col( source ) * it.value();
            }
            matrix.prune( 0, 0 );
        }

        // make the max row in idx-th column to be 1.
        void normalize( Index idx ) {
            SparseMatrix::ReverseInnerIterator it( matrix, idx );
            if ( it ) {
                matrix.col( idx ) /= it.value();
            }
        }

        void prune() {
            matrix.prune( 0, 0 );
        }

        // return set of nonzero indices at col
        std::set< Index > set_of_indices( const Index idx ) {
            std::set< Index > temp_set;
            matrix.prune( 0, 0 );
            for ( gyoza::ZpSparseMatrix::InnerIterator it( matrix, idx ); it; ++it ) {
                temp_set.insert( it.row());
            }
            return temp_set;
        }

        void set_matrix_value( const Index row, const Index col, const Core::Zp& val ) {
            matrix.coeffRef( row, col ) = val;
        }

        void print_matrix() {
            std::cout << matrix << std::endl;
        }

        // return idx-th set of vertices
        // Remark: call immediately after generating BoundaryMatrix
        IndexColumn get_vertices_of_simplex( const Index idx, const IndexMatrix& vertex_matrix ) {
            if ( this->is_zero( idx )) {
                IndexColumn atom;
                atom.push_back( idx );
                return atom;
            } else {
                auto subsimplex_set = this->set_of_indices( idx );
                auto itr = subsimplex_set.begin();
                IndexColumn vertices;
                if ( vertex_matrix.is_empty( *itr )) {
                    vertices = this->get_vertices_of_simplex( *itr, vertex_matrix );
                } else {
                    vertex_matrix.get_col( *itr, vertices );
                }
                ++itr;
                IndexColumn temp_vertices;
                if ( vertex_matrix.is_empty( *itr )) {
                    temp_vertices = this->get_vertices_of_simplex( *itr, vertex_matrix );
                } else {
                    vertex_matrix.get_col( *itr, temp_vertices );
                }
                vertices.insert( vertices.end(), temp_vertices.begin(), temp_vertices.end());
                std::sort( vertices.begin(), vertices.end());
                vertices.erase( std::unique( vertices.begin(), vertices.end()), vertices.end());
                return vertices;
            }
        }
    };

    class BoundaryMatrixDom : public BoundaryMatrix {
    protected:
        std::vector< Index > block;

    public:
        void resize_cols( Index cols ) {
            block.resize( cols );
            dims.resize( cols );
            matrix.resize( cols, cols );
        }

        void set_all_block( std::vector< Index > block_index ) {
            block = block_index;
        }

        void generate_new_bdmat_by_sorted_index( const IndexColumn& sorted_index,
                                                 const BoundaryMatrix& original_bdmat ) {
            const Index size = original_bdmat.cols();
            this->resize_cols( size );
            for ( Index original_idx = 0; original_idx != size; ++original_idx ) {
                //Index in dom filtration
                Index dom_idx = sorted_index.at( original_idx );

                //new column of the boundary matrix
                SparseColumn old_col, new_col;
                original_bdmat.get_col( original_idx, old_col );
                new_col.resize( old_col.size());
                for ( SparseColumn::InnerIterator it( old_col ); it; ++it ) {
                    new_col.coeffRef( sorted_index.at( it.row())) = it.value();
                }
                this->set_col( dom_idx, new_col );

                //copy dims
                this->set_dim( dom_idx, original_bdmat.get_dim( original_idx ));
            }
        }

        // sorted_domain_index = Index map of inclusion_map
        void set_boundary_matrix_dom( const IndexColumn& index_at_codomain, const BoundaryMatrix& original_bdmat,
                                      const Radii& radii_vector, const Radii& radii_vector_codomain,
                                      SparseMatrix& inclusion_map,
                                      IndexColumn& sorted_domain_index,
                                      CriticalIndexRadiusPairs& critical_index_and_radii ) {
            const Index size = original_bdmat.cols();
            this->resize_cols( size );
            std::vector< bool > flags_already_get( size, false ); // a vector f_a_g = (false, false, ..., false)

            IndexColumn sorted_index;// inverse map of sorted_domain_index

            sorted_domain_index.reserve( size );

            std::vector< Index > block_index;
            block_index.resize( size );

            //Index in idx-th dom, then Index <= idx && map(Index) <= filter_idx_codomain
            Radii unique_radii = radii_vector; // set of critical values (of domain and of codomain)
            unique_radii.insert( unique_radii.end(), radii_vector_codomain.begin(), radii_vector_codomain.end());
            sort_and_unique( unique_radii );
//            unique_radii.erase( std::unique( unique_radii.begin(), unique_radii.end()), unique_radii.end());
            int nr_critical_value = unique_radii.size();
            critical_index_and_radii.resize( nr_critical_value );
            Index critical_index_counter = -1;
            for ( Index i = 0; i < nr_critical_value; ++i ) {
                Index filter_idx_domain = max_index_under_value( radii_vector, unique_radii[i] );
                Index filter_idx_codomain = max_index_under_value( radii_vector_codomain, unique_radii[i] );
                for ( Index idx = 0; idx != filter_idx_domain + 1; ++idx ) {
                    if ( !flags_already_get.at( idx ) && index_at_codomain.at( idx ) != -1 &&
                         index_at_codomain.at( idx ) <= filter_idx_codomain ) {
                        sorted_domain_index.push_back( idx );
                        block_index[idx] = i;
                        flags_already_get.at( idx ) = true;
                        ++critical_index_counter;
                    }
                }
                critical_index_and_radii[i] = std::make_pair( critical_index_counter, unique_radii[i] );
            }
            // Fill -1 when image does not exist (never happen when cotractible)
            if ( critical_index_counter != size - 1 ) {
                for ( Index idx = 0; idx != size; ++idx ) {
                    if ( !flags_already_get.at( idx )) {
                        sorted_domain_index.push_back( idx );
                        block_index[idx] = -1;
                        flags_already_get.at( idx ) = true;
                        ++critical_index_counter;
                    }
                }
                critical_index_and_radii[nr_critical_value - 1] = std::make_pair( critical_index_counter, DBL_MAX );
            }

            // generate inclusion_map from sorted_domain_index
            inclusion_map.resize( size, size );
            for ( Index idx = 0; idx < size; ++idx ) {
                inclusion_map.coeffRef( sorted_domain_index.at( idx ), idx ) = 1;
            }
            inclusion_map.prune( 0, 0 );

            // generate sorted_index as the inverse of sorted_domain_index
            sorted_index.resize( sorted_domain_index.size());
            for ( Index idx = 0; idx < size; ++idx ) {
                sorted_index[sorted_domain_index.at( idx )] = idx;
            }

            // sort block_index according to sorted_index (this is equal to sort(block_index))
            std::vector< Index > temp_block_index;
            temp_block_index.resize( size );
            for ( Index idx = 0; idx < size; ++idx ) {
                temp_block_index[sorted_index[idx]] = block_index[idx];
            }
            std::swap( block_index, temp_block_index );

            if ( sorted_index.size() != size ) {
                std::cerr << "Function set_boundary_matrix_dom ERROR: sorted_index" << std::endl;
                throw;
            }

            if ( block_index.size() != size ) {
                std::cerr << "Function set_boundary_matrix_dom ERROR: block_index" << std::endl;
                throw;
            }

            this->generate_new_bdmat_by_sorted_index( sorted_index, original_bdmat );
            this->set_all_block( block_index );
        }
    };
}

#endif //PDSM_PERSISTENT_HOMOLOGY_H
