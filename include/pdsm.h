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

#ifndef PDSM_PDSM_H
#define PDSM_PDSM_H

#include "common_definitions_pdsm.h"
#include "persistent_homology.h"

namespace Pdsm {
    // Loads the map from given file
    // Format: i-th line represents an Index j. (i->j) of sampling data.
    bool load_map_vector( const std::string& filename, IndexColumn& map_vector ) {
        map_vector.clear();

        std::string cur_line;
        std::ifstream input_stream( filename.c_str());
        if ( input_stream.fail())
            return false;

        IndexColumn temp_col;
        while ( getline( input_stream, cur_line )) {
            if ( !cur_line.empty()) {
                std::stringstream ss( cur_line );
                int64_t temp_index;
                ss >> temp_index;
                map_vector.push_back((Index) temp_index );
            }
        }

        input_stream.close();
        return true;
    }

    void print_column( const IndexColumn& col ) {
        for ( const auto& idx : col ) {
            std::cout << idx << std::endl;
        }
    }

    void sort( IndexColumn& col ) {
        std::sort( col.begin(), col.end());
    }

    int sgn( const IndexColumn& col1, const IndexColumn& col2 ) {
        size_t size = col1.size();
        if ( size != col2.size()) {
            std::cerr << "Size Error: gyoza::sgn function." << std::endl;
            throw;
        }
        int counter = 0;
        for ( int i = 0; i < size; ++i ) {
            if ( col1[i] > col2[i] ) {
                ++counter;
            }
        }
        if ( counter % 2 == 0 ) { // even
            return 1;
        } else if ( counter % 2 == 1 ) { // odd
            return -1;
        } else {
            std::cerr << "Counter Value Error: gyoza::sgn function." << std::endl;
            throw;
        }
    }

    // generate index at codomain
    // Remark: call immediately after generating BoundaryMatrix
    void
    get_index_at_codomain( const IndexMatrix& vertex_matrix_of_codomain, const IndexMatrix& simplicial_map_vertex,
                           const boost::unordered_map< Pdsm::IndexColumn, Pdsm::Index >& sorted_vertices_codomain_to_index,
                           IndexColumn& index_at_codomain ) {
        const Index size = simplicial_map_vertex.cols(); // to get # of indices in domain
        index_at_codomain.resize( size );
        boost::progress_display show_progress( size );
        for ( Index idx = 0; idx != size; ++idx ) {
            ++show_progress;
            IndexColumn temp_col;
            simplicial_map_vertex.get_col( idx, temp_col );
            std::sort( temp_col.begin(), temp_col.end());
            temp_col.erase( std::unique( temp_col.begin(), temp_col.end()), temp_col.end());
            index_at_codomain[idx] = sorted_vertices_codomain_to_index.at( temp_col );
        }
    }

    bool non_degenerate( const IndexColumn& simplex ) {
        IndexColumn temp_simplex = simplex;

        std::sort( temp_simplex.begin(), temp_simplex.end());
        temp_simplex.erase( std::unique( temp_simplex.begin(), temp_simplex.end()), temp_simplex.end());

        return (temp_simplex.size() == simplex.size());
    }

    // From simplicial_map_vertex, make homological_simplicial_map
    // Ignore degenerate image.
    // Remark: call immediately after generating BoundaryMatrix
    void
    get_homological_simplicial_map( const IndexMatrix& vertex_matrix_codomain,
                                    const IndexMatrix& simplicial_map_vertex,
                                    const boost::unordered_map< Pdsm::IndexColumn, Pdsm::Index >& sorted_vertices_codomain_to_index,
                                    SparseMatrix& homological_simplicial_map ) {
        const Index size_domain = simplicial_map_vertex.cols(); // to get # of simplices in domain
        const Index size_codomain = vertex_matrix_codomain.cols();
        homological_simplicial_map.resize( size_codomain, size_domain );
        boost::progress_display show_progress( size_domain );
        for ( Index idx = 0; idx != size_domain; ++idx ) {
            ++show_progress;
            IndexColumn simplex;
            simplicial_map_vertex.get_col( idx, simplex );
            if ( non_degenerate( simplex )) {
                auto sorted_simplex = simplex;
                std::sort( sorted_simplex.begin(), sorted_simplex.end());
                Index idx_in_codomain = sorted_vertices_codomain_to_index.at( sorted_simplex );
                IndexColumn codomain_simplex;
                vertex_matrix_codomain.get_col( idx_in_codomain, codomain_simplex );
                Core::Zp coeff = sgn( simplex, codomain_simplex );

                homological_simplicial_map.coeffRef( idx_in_codomain, idx ) = coeff;
            }
        }
        homological_simplicial_map.prune( 0, 0 );
    }


    // make right projection q
    void
    get_right_projection( const SparseMatrix& left_projection,
                          const SparseMatrix& homological_simplicial_map,
                          SparseMatrix& right_projection ) {
        right_projection = homological_simplicial_map * left_projection;
    }


    class Reduced {
    private:
        BoundaryMatrix reduced_bdmat;

    protected:
        std::map< Index, Index > birth_death_pair;

    public:
        Index get_death_index_of_birth( Index birth ) const {
            return birth_death_pair.at( birth );
        }

        void get_col( Index idx, SparseColumn& col ) const {
            reduced_bdmat.get_col( idx, col );
        }

        // return generators = set of death indices s.t. idx \in (b,d) and dim dimensional.
        void get_generators_at( Index idx, Dimension dim, IndexColumn& generators ) const {
            generators.clear();

            for ( auto bd = birth_death_pair.begin(), bd_end = birth_death_pair.end(); bd != bd_end; ++bd ) {
                const Index birth = bd->first;
                const Index death = bd->second;
                if ( birth <= idx && idx < death && reduced_bdmat.get_dim( birth ) == dim ) {
                    generators.push_back( death );
                }
            }

            sort( generators );
        }

        void operator()( const BoundaryMatrix& bdmat ) {
            reduced_bdmat = bdmat;

            const Index cols = reduced_bdmat.cols();
            std::vector< Index > lowest( cols, -1 ); // a vector L = (-1, -1, ..., -1)

            for ( Index cur_col = 0; cur_col < cols; cur_col++ ) {
                Index max_idx = reduced_bdmat.get_max_index( cur_col );
                while ( max_idx != -1 && lowest[max_idx] != -1 ) {
                    reduced_bdmat.reduce( lowest[max_idx], cur_col ); // MATRIX ADDITION HERE!
                    max_idx = reduced_bdmat.get_max_index( cur_col );
                }
                if ( max_idx != -1 ) {
                    lowest[max_idx] = cur_col;
                    birth_death_pair.insert( std::make_pair( max_idx, cur_col ));
                    reduced_bdmat.normalize( cur_col ); // lowest non-zero value should be 1.
                }
                reduced_bdmat.prune();
            }
        }
    };

    class ReducedGraph : public Reduced {
    private:
        BoundaryMatrixGraph reduced_graph_bdmat;

    public:
        void get_col( Index idx, SparseColumn& col ) const {
            reduced_graph_bdmat.get_col( idx, col );
        }

        // return generators = set of death indices s.t. idx \in (b,d) and dim dimensional.
        void get_generators_at( Index idx, Dimension dim, IndexColumn& generators ) const {
            generators.clear();

            for ( const auto& bd : birth_death_pair ) {
                Index birth = bd.first; //birth at Index number.
                Index death = bd.second; //death at Index number.
                if ( birth <= idx && idx < death && reduced_graph_bdmat.get_dim( bd.first ) == dim ) {
                    generators.push_back( bd.second );
                }
            }

            sort( generators );
        }

        void operator()( const BoundaryMatrixGraph& bdmat_graph ) {
            reduced_graph_bdmat = bdmat_graph;

            const Index cols = reduced_graph_bdmat.cols();
            std::vector< Index > lowest( cols, -1 ); // a vector L = (-1, -1, ..., -1)

            for ( Index cur_col = 0; cur_col < cols; cur_col++ ) {
                Index max_idx = reduced_graph_bdmat.get_max_index( cur_col );
                while ( max_idx != -1 && lowest[max_idx] != -1 ) {
                    reduced_graph_bdmat.reduce( lowest[max_idx], cur_col ); // MATRIX ADDITION HERE!
                    max_idx = reduced_graph_bdmat.get_max_index( cur_col );
                }
                if ( max_idx != -1 ) {
                    lowest[max_idx] = cur_col;
                    birth_death_pair.insert( std::make_pair( max_idx, cur_col ));
                    reduced_graph_bdmat.normalize( cur_col ); // lowest non-zero value should be 1.
                }
                reduced_graph_bdmat.prune();
            }
        }
    };

    // STEP2
    // make homology induced maps as matrices

    // mapping cycles (non-zero SparseColumn in reduced_graph), and decompose by the basis of the target space
    void make_cycle_decomposition( const SparseColumn& cycle, const SparseMatrix& simplicial_map,
                                   const Reduced& reduce_target,
                                   SparseColumn& cycle_decomposition_of_image ) {
        // image_of_cycle:
        SparseColumn image_of_cycle = simplicial_map * cycle;
        image_of_cycle.prune( 0, 0 );

        while ( image_of_cycle.nonZeros() != 0 ) {
            SparseColumn::ReverseInnerIterator it( image_of_cycle );
            const Index max_idx = it.row();
            Index source_col_idx = reduce_target.get_death_index_of_birth( max_idx );

            SparseColumn basis;
            reduce_target.get_col( source_col_idx, basis );

            cycle_decomposition_of_image.coeffRef( source_col_idx ) = it.value();
            image_of_cycle -= basis * it.value();

            image_of_cycle.prune( 0, 0 );
        }
    }

    // compute homology induced map matrix (cycle_matrix)
    void make_cycle_matrix( const ReducedGraph& reduce_graph, const SparseMatrix& simplicial_map,
                            const Reduced& reduce_target,
                            SparseMatrix& cycle_matrix ) {
        cycle_matrix.resize( simplicial_map.rows(), simplicial_map.cols());
        auto col_size = simplicial_map.cols();
        for ( Index idx = 0; idx < col_size; ++idx ) {
            SparseColumn cycle;
            reduce_graph.get_col( idx, cycle );

            SparseColumn cycle_decomposition_of_image;
            cycle_decomposition_of_image.resize( simplicial_map.cols());
            cycle.prune( 0, 0 );
            if ( cycle.nonZeros() != 0 ) {
                make_cycle_decomposition( cycle, simplicial_map, reduce_target, cycle_decomposition_of_image );
            }

            cycle_matrix.col( idx ) = cycle_decomposition_of_image;
        }
        cycle_matrix.prune( 0, 0 );
    }

    // get idx-th matrix from homological_map(left projection p OR right projection q)
    void make_cycle_matrix_at( Index idx, Radius rad, Dimension dim, const ReducedGraph& reduce_graph,
                               const Reduced& reduce_target, const Radii& radii_vector_target,
                               const SparseMatrix& cycle_matrix, DenseMatrix& cycle_matrix_at_idx ) {
        // generators of HG_{idx}
        IndexColumn generators_graph_at_idx;
        reduce_graph.get_generators_at( idx, dim, generators_graph_at_idx );
        gyoza::Index cols = generators_graph_at_idx.size();

        // generators of HC_{maxidx <= rad} or HD_{maxidx <= rad}
        Index idx_target = max_index_under_value( radii_vector_target, rad );
        IndexColumn generators_at_idx;
        reduce_target.get_generators_at( idx_target, dim, generators_at_idx );
        gyoza::Index rows = generators_at_idx.size();

        cycle_matrix_at_idx = DenseMatrix::Zero( rows, cols );

        // in Matlab, this is cycle_matrix[generators_at_idx, generators_graph_at_idx]
        for ( gyoza::Index col_idx = 0; col_idx < cols; ++col_idx ) {
            SparseColumn temp_col = cycle_matrix.col( generators_graph_at_idx[col_idx] );
            for ( gyoza::Index row_idx = 0; row_idx < rows; ++row_idx ) {
                cycle_matrix_at_idx.coeffRef( row_idx, col_idx ) = temp_col.coeff( generators_at_idx[row_idx] );
            }
        }
    }

    // maps in HG filtration (denoted as filter_matrix）
    void make_filter_matrix_between( Index idx_source, Index idx_target, Dimension dim,
                                     const ReducedGraph& reduce_graph,
                                     DenseMatrix& filter_matrix_between_idx ) {
        // generators of HG_{idx_source}
        IndexColumn generators_source;
        reduce_graph.get_generators_at( idx_source, dim, generators_source );
        gyoza::Index cols = generators_source.size();

        // generators of HG_{idx_target}
        IndexColumn generators_target;
        reduce_graph.get_generators_at( idx_target, dim, generators_target );
        gyoza::Index rows = generators_target.size();

        filter_matrix_between_idx = DenseMatrix::Zero( rows, cols );

        for ( Index jj = 0; jj < cols; ++jj ) {
            for ( Index ii = 0; ii < rows; ++ii ) {
                if ( generators_source.at( jj ) == generators_target.at( ii )) {
                    filter_matrix_between_idx( ii, jj ) = 1;
                }
            }
        }
    }


    // return: dim of I[1,3] computed from p_* and q_*
    gyoza::Index get_change_of_basis( const DenseMatrix& induced_left_matrix, const DenseMatrix& induced_right_matrix,
                                      gyoza::algorithms::FullPivotSNF< DenseMatrix >& snf_base_change ) {
        auto full_cols = induced_left_matrix.cols();
        DenseMatrix graph_base_change = DenseMatrix::Identity( full_cols, full_cols );

        gyoza::algorithms::FullPivotSNF< DenseMatrix > snf_inc( induced_left_matrix, true );

        // base change
        graph_base_change = graph_base_change * snf_inc.get_Q();

        auto rank_inc = snf_inc.rank();

        // Note: q_matrix is the induced right matrix Q_*
        DenseMatrix sub_left_q_matrix = (induced_right_matrix * snf_inc.get_Q()).leftCols( rank_inc );
        DenseMatrix sub_right_q_matrix = (induced_right_matrix * snf_inc.get_Q()).rightCols( full_cols - rank_inc );
        gyoza::algorithms::FullPivotSNF< DenseMatrix > snf_right_q( sub_right_q_matrix, true );
        sub_left_q_matrix = snf_right_q.get_P() * sub_left_q_matrix;

        // base change
        DenseMatrix temp_base_change = DenseMatrix::Identity( full_cols, full_cols );
        temp_base_change.bottomRightCorner( snf_right_q.get_Q().rows(),
                                            snf_right_q.get_Q().cols()) = snf_right_q.get_Q();
        graph_base_change = graph_base_change * temp_base_change;

        auto rank_right_q = snf_right_q.rank();

        DenseMatrix sub_left_top_q_matrix = sub_left_q_matrix.topRows( rank_right_q );
        DenseMatrix sub_left_bottom_q_matrix = sub_left_q_matrix.bottomRows(
                induced_right_matrix.rows() - rank_right_q );

        // base change (zeroing out sub_left_top_q_matrix by snf of sub_right_q)
        temp_base_change = DenseMatrix::Identity( full_cols, full_cols );
        temp_base_change.block( sub_left_top_q_matrix.cols(), 0, sub_left_top_q_matrix.rows(),
                                sub_left_top_q_matrix.cols()) = (-1) * sub_left_top_q_matrix;
        graph_base_change = graph_base_change * temp_base_change;

        gyoza::algorithms::FullPivotSNF< DenseMatrix > snf_left_bottom_q( sub_left_bottom_q_matrix, true );

        // base change
        temp_base_change = DenseMatrix::Identity( full_cols, full_cols );
        temp_base_change.topLeftCorner( snf_left_bottom_q.get_Q().rows(),
                                        snf_left_bottom_q.get_Q().cols()) = snf_left_bottom_q.get_Q();
        graph_base_change = graph_base_change * temp_base_change;

        snf_base_change.compute( graph_base_change );//Using snf to get inverse

        return snf_left_bottom_q.rank(); //dim of I[1,3]
    }

    // make the main (last) persistence module
    void
    get_all_change_of_basis( const Dimension dim, const ReducedGraph& reduce_graph,
                             const Reduced& reduce_domain,
                             const Reduced& reduce_codomain,
                             const Radii& radii_vector_domain,
                             const Radii& radii_vector_codomain,
                             const SparseMatrix& induced_left_matrix, const SparseMatrix& induced_right_matrix,
                             const CriticalIndexRadiusPairs& critical_index_radius_pairs,
                             std::vector< DenseMatrix >& main_persistence_module,
                             std::vector< DenseMatrix >& change_basis_generators ) {
        const Index nr_critical_values = critical_index_radius_pairs.size();
        change_basis_generators.resize( nr_critical_values );
        main_persistence_module.resize( nr_critical_values - 1 );
        std::vector< gyoza::Index > module_dims;
        module_dims.resize( nr_critical_values );
        for ( Index block = 0; block < nr_critical_values; ++block ) {
            const Index critical_idx = critical_index_radius_pairs.at( block ).first; // max index in block
            const Radius critical_radius = critical_index_radius_pairs.at( block ).second;
            // p_*_{critical_radius}
            DenseMatrix induced_left_matrix_at_block;
            make_cycle_matrix_at( critical_idx, critical_radius, dim,
                                  reduce_graph, reduce_domain, radii_vector_domain,
                                  induced_left_matrix,
                                  induced_left_matrix_at_block );
            // q_*_{critical_radius}
            DenseMatrix induced_right_matrix_at_block;
            make_cycle_matrix_at( critical_idx, critical_radius, dim,
                                  reduce_graph, reduce_codomain, radii_vector_codomain,
                                  induced_right_matrix,
                                  induced_right_matrix_at_block );

            // H(inclusion of graph filt) from critical_idx to next critical_idx
            if ( block < nr_critical_values - 1 ) {
                Index idx_source = critical_index_radius_pairs[block].first;
                Index idx_target = critical_index_radius_pairs[block + 1].first;
                make_filter_matrix_between( idx_source, idx_target, dim, reduce_graph,
                                            main_persistence_module[block] );
            }

            // compute basis of I[1,3] in HG_idx
            gyoza::algorithms::FullPivotSNF< DenseMatrix > snf_base_change;
            module_dims[block] = get_change_of_basis( induced_left_matrix_at_block, induced_right_matrix_at_block,
                                                      snf_base_change );


            DenseMatrix base_change = snf_base_change.get_mat();
            if ( !snf_base_change.get_inverse()) {
                std::cerr << "Caution: snf_base_change is not full rank!" << std::endl;
                std::cerr << base_change << std::endl;
                std::cerr << block << std::endl;
                std::cerr << module_dims.at( block ) << std::endl;
                std::cerr << induced_left_matrix_at_block << std::endl << std::endl;
                std::cerr << induced_right_matrix_at_block << std::endl;
                throw;
            }
            DenseMatrix base_change_inv = *snf_base_change.get_inverse();

            // change the basis of and truncate morphism in persistence module
            // as target
            if ( block > 0 ) {
                main_persistence_module[block - 1] = (base_change_inv *
                                                      main_persistence_module[block - 1]).topRows(
                        module_dims[block] );
            }
            // as source
            if ( block < nr_critical_values - 1 ) {
                main_persistence_module[block] = (main_persistence_module[block] * base_change).leftCols(
                        module_dims[block] );
            }

            change_basis_generators[block] = base_change *
                                             DenseMatrix::Identity( base_change.cols(), module_dims[block] );
        }

    }


    // transform column echelon form
    void transform_col_echelon( DenseMatrix& mat, DenseMatrix& col_base_change_inv,
                                DenseMatrix& change_basis_gen_here,
                                std::vector< gyoza::Index >& pivots ) {
        const gyoza::Index col_size = mat.cols();
        const gyoza::Index row_size = mat.rows();

        DenseMatrix col_base_change = DenseMatrix::Identity( col_size, col_size );

        gyoza::Index pivot_counter = 0;

        for ( gyoza::Index piv_row = 0; piv_row < row_size; ++piv_row ) {
            for ( gyoza::Index piv_col = pivot_counter; piv_col < col_size; ++piv_col ) {
                if ( mat( piv_row, piv_col ) != 0 ) { // pivot found
                    // normalize
                    mat.col( piv_col ) /= mat( piv_row, piv_col );
                    col_base_change( piv_col ) /= mat( piv_row, piv_col );
                    // eliminate nonzero values rightward of the pivot using pivot
                    for ( gyoza::Index cur_col = piv_col + 1; cur_col < col_size; ++cur_col ) {
                        if ( mat( piv_row, cur_col ) != 0 ) {
                            mat.col( cur_col ) -= mat( piv_row, cur_col ) * mat.col( piv_col );
                            col_base_change.col( cur_col ) -=
                                    mat( piv_row, cur_col ) * col_base_change.col( piv_col );
                        }
                    }
                    mat.col( pivot_counter ).swap( mat.col( piv_col ));
                    col_base_change.col( pivot_counter ).swap( col_base_change.col( piv_col ));
                    pivots.push_back( piv_row );
                    ++pivot_counter;
                    break;
                }
            }
        }

//        // compute inverse change of basis
        gyoza::algorithms::FullPivotSNF< DenseMatrix > snf_base_change( col_base_change, true );
        col_base_change_inv = *snf_base_change.get_inverse();

        // store the change of basis to calculate generators by backtrace
        change_basis_gen_here = change_basis_gen_here * col_base_change;
    }

    // eliminate elements lower than each pivot in mat
    // preserve matrices which give elementary operations
    void eliminate_lower_pivot( DenseMatrix& mat, const std::vector< gyoza::Index >& pivots,
                                DenseMatrix& change_basis_generators_next_block ) {
        size_t nr_pivots = pivots.size();
        long nr_rows = mat.rows();

        DenseMatrix row_base_change = DenseMatrix::Identity( nr_rows, nr_rows );

        for ( auto piv_col = 0; piv_col < nr_pivots; ++piv_col ) {
            gyoza::Index piv_row = pivots.at( piv_col );
            for ( auto cur_row = piv_row + 1; cur_row < nr_rows; ++cur_row ) {
                if ( mat( cur_row, piv_col ) != 0 ) {
                    mat.row( cur_row ) -= mat( cur_row, piv_col ) * mat.row( piv_row );
                    row_base_change.row( cur_row ) -= mat( cur_row, piv_col ) * row_base_change.row( piv_row );
                }
            }
        }

        // compute inverse change of basis
        gyoza::algorithms::FullPivotSNF< DenseMatrix > snf_base_change( row_base_change, true );
        DenseMatrix row_base_change_inv = *snf_base_change.get_inverse();
        change_basis_generators_next_block = change_basis_generators_next_block * row_base_change_inv;
    }

    // class preserving birth-death pairs
    class BirthDeathPairs {
    private:
        std::vector< std::pair< Index, Index >> pairs;

    public:
        Index get_num_pairs() const {
            return pairs.size();
        }

        void append_pair( Index birth, Index death ) {
            pairs.push_back( std::make_pair( birth, death ));
        }

        std::pair< Index, Index > get_pair_at( Index idx ) const {
            return pairs[idx];
        }
    };

    // make persistence diagram (ref. sec 3.4 Tower Bases in the original paper)
    void
    make_persistence_diagram( const Dimension dim,
                              std::vector< DenseMatrix >& main_persistence_module,
                              const CriticalIndexRadiusPairs& critical_index_radius_pairs,
                              std::vector< DenseMatrix >& change_basis_generators,
                              const ReducedGraph& reduce_graph,
                              BirthDeathPairs& bd_pairs,
                              std::vector< gyoza::Index >& piv_col_of_bd_pairs,
                              const IndexColumn& sorted_graph_index,
                              const IndexMatrix& vertex_matrix,
                              const std::string& outfile_directory ) {
        Index nr_blocks = main_persistence_module.size() + 1; //size = # of maps
        std::vector< std::vector< gyoza::Index > > vec_pivots;
        vec_pivots.resize( main_persistence_module.size());

        // to column echelon form
        std::cout << "\nComputing column echelon forms:";
        boost::progress_display show_progress_col_echelon( nr_blocks - 1 );
        for ( Index cur_block = nr_blocks - 2; cur_block > -1; --cur_block ) {
            ++show_progress_col_echelon;
            DenseMatrix col_base_change_inv;
            transform_col_echelon( main_persistence_module[cur_block], col_base_change_inv,
                                   change_basis_generators[cur_block],
                                   vec_pivots[cur_block] );

            //side effects to the lower next filter
            if ( cur_block != 0 ) {
                main_persistence_module[cur_block - 1] =
                        col_base_change_inv * main_persistence_module[cur_block - 1];
            }
        }

        // eliminate elements lower than pivots
        std::cout << "\nComputing row echelon forms:";
        boost::progress_display show_progress_row_echelon( nr_blocks - 1 );
        for ( Index cur_block = 0; cur_block < nr_blocks - 1; ++cur_block ) {
            ++show_progress_row_echelon;
//            DenseMatrix row_base_change_inv;
            eliminate_lower_pivot( main_persistence_module[cur_block], vec_pivots[cur_block],
                                   change_basis_generators[cur_block + 1] );
        }

        // get birth-death pairs using pivots
        // Remark: Index here is block number
        std::vector< std::pair< Index, Index > > birth_val_pairs;
        std::cout << "\nComputing birth death pairs:";
        std::vector< gyoza::Index > temp_piv_col; // record which Index at birth
        boost::progress_display show_progress( nr_blocks );
        // initialize
        auto nr_0th_pivots = vec_pivots.at( 0 ).size();
        for ( int idx = 0; idx < nr_0th_pivots; ++idx ) {
            birth_val_pairs.emplace_back( 0, idx );
            temp_piv_col.push_back( idx );
        }
        // basis of Ker gives birth-death pair of [0,1].
        auto dim_at0 = main_persistence_module[0].cols();
        for ( int idx = nr_0th_pivots; idx < dim_at0; ++idx ) {
            bd_pairs.append_pair( 0, 1 );
            piv_col_of_bd_pairs.push_back( idx );
        }
        ++show_progress;
        // after 1st filter
        for ( auto cur_block = 1; cur_block < nr_blocks - 1; ++cur_block ) {
            ++show_progress;
            // Last pivots give infinity-death // TODO: implement inf value (never happen when contractible)
            if ( cur_block == nr_blocks - 2 ) {
                auto nr_pairs = birth_val_pairs.size();
                auto nr_pivots = vec_pivots.at( cur_block ).size();
                std::vector< Index > kernel_indices_target_from_former_pivots;
                for ( int idx = 0; idx < nr_pairs; ++idx ) {
                    auto pair = birth_val_pairs[idx];
                    // We get death values ... (Ker from cur_block to cur_block + 1)
                    if ( vec_pivots.at( cur_block - 1 ).at( pair.second ) >= nr_pivots ) {
                        bd_pairs.append_pair( pair.first, cur_block + 1 );
                        piv_col_of_bd_pairs.push_back( temp_piv_col[idx] );
                        kernel_indices_target_from_former_pivots.push_back(
                                vec_pivots.at( cur_block - 1 ).at( pair.second ));
                    } else { // never happen when contractible
                        bd_pairs.append_pair( pair.first, cur_block + 2 );
                        piv_col_of_bd_pairs.push_back( temp_piv_col[idx] );
                    }
                    // intervals [cur_block, cur_block+1]
                    auto dim_source = main_persistence_module[cur_block].cols();
                    auto rank = vec_pivots[cur_block].size();
                    for ( int idx = rank; idx < dim_source; ++idx ) {
                        auto itr = std::find( kernel_indices_target_from_former_pivots.begin(),
                                              kernel_indices_target_from_former_pivots.end(), idx );
                        if ( itr == kernel_indices_target_from_former_pivots.end()) {
                            bd_pairs.append_pair( cur_block, cur_block + 1 );
                            piv_col_of_bd_pairs.push_back( idx );
                        }
                    }
                }
            } else {
                std::vector< std::pair< Index, Index > > temp_pairs;
                std::vector< gyoza::Index > temp_temp_piv_col;
                auto nr_pivots = vec_pivots.at( cur_block ).size();
                std::vector< Index > used_indices;
                std::vector< Index > kernel_indices_target_from_former_pivots;
                auto nr_pairs = birth_val_pairs.size();
                for ( int idx = 0; idx < nr_pairs; ++idx ) {
                    auto pair = birth_val_pairs[idx];

                    // if greater than nr_pivots, the next filter is death value (Ker from cur_block to cur_block + 1)
                    if ( vec_pivots.at( cur_block - 1 ).at( pair.second ) >= nr_pivots ) {
                        bd_pairs.append_pair( pair.first, cur_block + 1 );
                        piv_col_of_bd_pairs.push_back( temp_piv_col[idx] );
                        kernel_indices_target_from_former_pivots.push_back(
                                vec_pivots.at( cur_block - 1 ).at( pair.second ));
                    } else {
                        temp_pairs.emplace_back( pair.first, vec_pivots.at( cur_block - 1 ).at( pair.second ));
                        temp_temp_piv_col.push_back( temp_piv_col[idx] );
                        used_indices.push_back( vec_pivots.at( cur_block - 1 ).at( pair.second ));
                    }
                }
                // How to get [cur_block, cur_block+1]...
                // #Ker = #[*, cur_block+1] = dim_source - rank
                // [cur_block, cur_block+1] = #Ker - #(generators out of range at former pivots)
                auto dim_source = main_persistence_module[cur_block].cols();
                auto rank = vec_pivots[cur_block].size();
                for ( int idx = rank; idx < dim_source; ++idx ) {
                    auto itr = std::find( kernel_indices_target_from_former_pivots.begin(),
                                          kernel_indices_target_from_former_pivots.end(), idx );
                    if ( itr == kernel_indices_target_from_former_pivots.end()) {
                        bd_pairs.append_pair( cur_block, cur_block + 1 );
                        piv_col_of_bd_pairs.push_back( idx );
                    }
                }
                for ( auto idx = 0; idx < nr_pivots; ++idx ) {
                    auto itr = std::find( used_indices.begin(), used_indices.end(), idx );
                    if ( itr == used_indices.end()) {
                        temp_pairs.emplace_back( cur_block, idx );
                        temp_temp_piv_col.push_back( idx );
                    }
                }
                birth_val_pairs.swap( temp_pairs );
                temp_piv_col.swap( temp_temp_piv_col );
            }
        }


        //print bd pairs
        auto nr_bd = bd_pairs.get_num_pairs();
        std::cout << "\nThe birth death pairs are:\n";
        for ( auto i = 0; i < nr_bd; ++i ) {
            auto pair = bd_pairs.get_pair_at( i );
            auto birth_radius = critical_index_radius_pairs.at( pair.first ).second;
            auto death_radius = critical_index_radius_pairs.at( pair.second ).second;
            std::ofstream outputfile_bdpair( outfile_directory + "/bd" + std::to_string( i ) + ".txt" );
            std::cout << "("
                      << birth_radius
                      << ", "
                      << death_radius
                      << ") Index:("
                      << critical_index_radius_pairs.at( pair.first ).first
                      << ", "
                      << critical_index_radius_pairs.at( pair.second ).first
                      << ") filter num:("
                      << pair.first
                      << ", "
                      << pair.second
                      << ")\n";
            outputfile_bdpair << birth_radius << " " << death_radius << "\n";
            outputfile_bdpair.close();

            std::cout << "generator:\n";

            //make canonical basis at birth
            Index birth_block_idx = pair.first;

            // canonical basis at pivot is a generator, so transform this by inverses operations.
            DenseMatrix can_basis = DenseMatrix::Zero( main_persistence_module[birth_block_idx].cols(), 1 );
            can_basis( piv_col_of_bd_pairs[i], 0 ) = 1;

            std::cout << can_basis << std::endl;
            std::cout << change_basis_generators[birth_block_idx] << std::endl;

            //transform it to the vector at filter matrix at b
            can_basis = change_basis_generators[birth_block_idx] * can_basis;

            //compute the generators of this cycle
            IndexColumn generators;
            reduce_graph.get_generators_at( critical_index_radius_pairs.at( birth_block_idx ).first, dim,
                                            generators );
            std::ofstream outputfile_generator( outfile_directory + "/generator" + std::to_string( i ) + ".txt" );
            for ( int j = 0; j < generators.size(); ++j ) {
                if ( can_basis( j, 0 ) != 0 ) {
                    std::cout << generators[j] << "\n";
                    std::cout << "coeff = " << can_basis( j, 0 ) << "\n";
                    SparseColumn vec_edges;
                    reduce_graph.get_col( generators[j], vec_edges );
                    for ( SparseColumn::InnerIterator it( vec_edges ); it; ++it ) {
                        auto edge_idx = it.row();
                        std::cout << "coeff = " << it.value() << ": ";
                        IndexColumn vec_vertices;
                        vertex_matrix.get_col( sorted_graph_index[edge_idx], vec_vertices );
                        std::cout << "(";
                        auto nr_vertices = vec_vertices.size();
                        for ( auto k = 0; k < nr_vertices; ++k ) {
                            auto vertex = vec_vertices[k];
                            if ( k != 0 ) {
                                std::cout << ", ";
                                outputfile_generator << " ";
                            }
                            std::cout << vertex;
                            outputfile_generator << vertex;
                        }
                        std::cout << ")\n";
                        outputfile_generator << "\n";
                    }
                }
            }
            std::cout << std::endl;
            outputfile_generator.close();
        }
    }
}

#endif //PDSM_PDSM_H
