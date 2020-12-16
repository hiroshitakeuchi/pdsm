//
// Created by Hiroshi Takeuchi on 2019/11/27.
//

/*  Copyright (C) 2019 Hiroshi Takeuchi

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

#include "include/pdsm.h"
#include "include/rips_torus.h"
#include <filesystem>

constexpr int16_t DIM = 1; // dimension of homology

int main( int argc, char *argv[] ) {
    Pdsm::Dimension dim = DIM;
    std::string infile_points_domain( argv[1] );
    std::string infile_points_codomain( argv[2] );
    std::string infile_map( argv[3] );
    std::string outfile_directory( argv[4] );

    auto time_filt = std::chrono::system_clock::now();
    // ########## STEP 0 making filtration ##########
    std::cout << "\n########## STEP 0 ##########" << std::endl;
    Pdsm::IndexColumn map_vector;
    Pdsm::load_map_vector( infile_map, map_vector );
    std::cout << "map data:" << std::endl;
    Pdsm::print_column( map_vector );

    Rips::PointContainer points_values_domain;
    Rips::IndexContainer to_new_index_domain;
    Pdsm::BoundaryMatrix bdmat_domain;
    Pdsm::Radii radii_vector_domain;
    std::vector <Rips::IndexContainer> complex_domain;
    Rips::create_rips_filtration( infile_points_domain, points_values_domain, to_new_index_domain, bdmat_domain,
                                  radii_vector_domain, complex_domain );
    bool flag_print = (bdmat_domain.cols() < 26);

    // boundary matrix of codomain
    Rips::PointContainer points_values_codomain;
    Rips::IndexContainer to_new_index_codomain;
    Pdsm::BoundaryMatrix bdmat_codomain;
    Pdsm::Radii radii_vector_codomain;
    std::vector <Rips::IndexContainer> complex_codomain;
    Rips::create_rips_filtration( infile_points_codomain, points_values_codomain, to_new_index_codomain, bdmat_codomain,
                                  radii_vector_codomain, complex_codomain );


    // ########## STEP0.5 create maps ##########
    std::cout << "\n########## STEP 0.5 ##########" << std::endl;
    std::cout << "The vertex matrix (orientations of simplices) of domain is" << std::endl;
    Pdsm::IndexMatrix vertex_matrix_domain;
    boost::unordered_map <Pdsm::IndexColumn, Pdsm::Index> sorted_vertices_domain_to_index;
    Rips::create_vertex_matrix( to_new_index_domain, complex_domain, vertex_matrix_domain,
                                sorted_vertices_domain_to_index );
    if ( flag_print ) {
        vertex_matrix_domain.print_matrix();
    }

    // create map of codomain
    Pdsm::IndexMatrix vertex_matrix_codomain;
    boost::unordered_map <Pdsm::IndexColumn, Pdsm::Index> sorted_vertices_codomain_to_index;
    Rips::create_vertex_matrix( to_new_index_codomain, complex_codomain, vertex_matrix_codomain,
                                sorted_vertices_codomain_to_index );


    std::cout << "\nThe simplicial map (as vertices) is:" << std::endl;
    Pdsm::IndexMatrix simplicial_map_vertex;
    vertex_matrix_domain.get_simplicial_map_vertex( map_vector, simplicial_map_vertex );
    if ( flag_print ) {
        simplicial_map_vertex.print_matrix();
    }

    std::cout << "\nThe Index at codomain is:" << std::endl;
    Pdsm::IndexColumn index_at_codomain;
    Pdsm::get_index_at_codomain( vertex_matrix_codomain, simplicial_map_vertex, sorted_vertices_codomain_to_index,
                                  index_at_codomain );
    if ( flag_print ) {
        Pdsm::print_column( index_at_codomain );
    }

    std::cout << "\nThe homological simplicial map is:" << std::endl;
    Pdsm::SparseMatrix homological_simplicial_map;
    Pdsm::get_homological_simplicial_map( vertex_matrix_codomain, simplicial_map_vertex,
                                           sorted_vertices_codomain_to_index, homological_simplicial_map );
    if ( flag_print ) {
        std::cout << homological_simplicial_map << std::endl;
    }

    std::cout << "\n### Generating graph boundary matrix ###" << std::endl;
    Pdsm::BoundaryMatrixGraph graph_bdmat;
    Pdsm::SparseMatrix left_projection;
    Pdsm::IndexColumn sorted_domain_index;
    Pdsm::CriticalIndexRadiusPairs critical_index_radius_pairs;
    graph_bdmat.set_boundary_matrix_graph( index_at_codomain, bdmat_domain, radii_vector_domain,
                                           radii_vector_codomain,
                                           left_projection, sorted_domain_index, critical_index_radius_pairs );
    std::cout << "\nThe graph boundary matrix is:" << std::endl;
    if ( flag_print ) {
        graph_bdmat.print_matrix();
    }
    std::cout << "\nThe left projection is:" << std::endl;
    if ( flag_print ) {
        std::cout << left_projection << std::endl;
    }

    std::cout << "\nThe right projection is:" << std::endl;
    Pdsm::SparseMatrix right_projection;
    Pdsm::get_right_projection( left_projection, homological_simplicial_map,
                                right_projection );
    if ( flag_print ) {
        std::cout << right_projection << std::endl;
    }

    auto time_start = std::chrono::system_clock::now();
    std::cout << "STEP0: "
              << std::chrono::duration_cast< std::chrono::milliseconds >( time_start - time_filt ).count() / 1000.
              << "s";
    // ########## STEP1 compute reduced matrix ##########
    std::cout << "\n########## STEP 1 ##########\nThe reduced boundary matrix is:" << std::endl;
    Pdsm::Reduced reduce;
    reduce( bdmat_domain );
    if ( flag_print ) {
        bdmat_domain.print_matrix();
    }

    // reduced matrix of codomain
    std::cout << "The reduced codomain boundary matrix is:" << std::endl;
    Pdsm::Reduced reduce_codomain;
    reduce_codomain( bdmat_codomain );

    std::cout << "The reduced graph boundary matrix is:" << std::endl;
    Pdsm::ReducedGraph reduce_graph;
    reduce_graph( graph_bdmat );
    if ( flag_print ) {
        graph_bdmat.print_matrix();
    }

    auto time_step1 = std::chrono::system_clock::now();
    std::cout << "STEP1: "
              << std::chrono::duration_cast< std::chrono::milliseconds >( time_step1 - time_start ).count() / 1000.
              << "s";
    // ########## STEP2 compute Matrix pairs ##########
    std::cout << "\n########## STEP 2 ##########\ncompute Matrix pairs" << std::endl;

    std::cout << "\nThe induced left matrix p_* is:" << std::endl;
    Pdsm::SparseMatrix induced_left_matrix;
    Pdsm::make_cycle_matrix( reduce_graph, left_projection, reduce, induced_left_matrix );
    if ( flag_print ) {
        std::cout << induced_left_matrix << std::endl;
    }

    std::cout << "\nThe induced right matrix q_* is:" << std::endl;
    Pdsm::SparseMatrix induced_right_matrix;
    Pdsm::make_cycle_matrix( reduce_graph, right_projection, reduce_codomain, induced_right_matrix );
    if ( flag_print ) {
        std::cout << induced_right_matrix << std::endl;
    }


    auto time_step2 = std::chrono::system_clock::now();
    std::cout << "STEP2: "
              << std::chrono::duration_cast< std::chrono::milliseconds >( time_step2 - time_step1 ).count() / 1000.
              << "s";
    // ########## STEP3 compute bases for I[1,3] ##########
    // compute the last persistence diagram
    std::cout << "\n########## STEP3 ##########\nComputing persistence module\n";
    std::vector <Pdsm::DenseMatrix> main_persistence_module;
    std::vector <Pdsm::DenseMatrix> change_basis_generators;
    Pdsm::get_all_change_of_basis( dim, reduce_graph,
                                   reduce, reduce_codomain,
                                   radii_vector_domain, radii_vector_codomain,
                                   induced_left_matrix, induced_right_matrix,
                                   critical_index_radius_pairs,
                                   main_persistence_module,
                                   change_basis_generators );


    std::cout << "\nComputing persistence diagram:\n";
    std::filesystem::create_directories( outfile_directory );
    Pdsm::BirthDeathPairs bd_pairs;
    std::vector <gyoza::Index> piv_col_of_bd_pairs;
    Pdsm::make_persistence_diagram( dim, main_persistence_module, critical_index_radius_pairs, change_basis_generators,
                                    reduce_graph,
                                    bd_pairs, piv_col_of_bd_pairs, sorted_domain_index, vertex_matrix_domain,
                                    outfile_directory );

    auto time_step3 = std::chrono::system_clock::now();
    std::cout << "STEP3: "
              << std::chrono::duration_cast< std::chrono::milliseconds >( time_step3 - time_step2 ).count() / 1000.
              << "s";

    return 0;
}