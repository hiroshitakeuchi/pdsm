//
// Created by Hiroshi Takeuchi on 2019/11/27.
// Difference from rips.h is only the distance function
// Support only 2-dim torus [0,1) * [0,1)
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

#ifndef PDSM_RIPS_TORUS_H
#define PDSM_RIPS_TORUS_H


#include "pdsm.h"

namespace Rips {
    typedef double DistanceType;
    typedef std::vector< DistanceType > DistanceContainer;
    typedef std::vector< double > Point;
    typedef std::vector< Point > PointContainer;
    typedef int64_t Index;
    typedef std::vector< Index > IndexContainer;
    typedef int Coeff;


    void read_points( const std::string& infilename, PointContainer& points ) {
        std::ifstream in( infilename.c_str());
        std::string line;
        while ( std::getline( in, line )) {
            std::stringstream linestream( line );
            double x;
            points.push_back( Point());
            while ( linestream >> x ) {
                points.back().push_back( x );
            }
        }
    }

    DistanceType l2distance( const Point& p1, const Point& p2 ) {
        DistanceType sum = 0;
        for ( size_t i = 0; i < p1.size(); ++i )
            sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);

        return sqrt( sum );
    }

    DistanceType distance( const Point& p1, const Point& p2 ) {
        // CAUTION: only considering p1, p2 \in [0,1) \times [0,1)
        // Method: copy 8 mirror points (right, rightup, up, upleft, left, leftdown, down, downright) and get the minimum
        DistanceType dist = l2distance( p1, p2 );
        DistanceType temp_dist;

        Point p1_right = p1;
        p1_right[0] += 1;
        temp_dist = l2distance( p1_right, p2 );
        dist = std::min( dist, temp_dist );

        Point p1_rightup = p1;
        p1_rightup[0] += 1;
        p1_rightup[1] += 1;
        temp_dist = l2distance( p1_rightup, p2 );
        dist = std::min( dist, temp_dist );

        Point p1_up = p1;
        p1_up[1] += 1;
        temp_dist = l2distance( p1_up, p2 );
        dist = std::min( dist, temp_dist );

        Point p1_upleft = p1;
        p1_upleft[0] -= 1;
        p1_upleft[1] += 1;
        temp_dist = l2distance( p1_upleft, p2 );
        dist = std::min( dist, temp_dist );

        Point p1_left = p1;
        p1_left[0] -= 1;
        temp_dist = l2distance( p1_left, p2 );
        dist = std::min( dist, temp_dist );

        Point p1_leftdown = p1;
        p1_leftdown[0] -= 1;
        p1_leftdown[1] -= 1;
        temp_dist = l2distance( p1_leftdown, p2 );
        dist = std::min( dist, temp_dist );

        Point p1_down = p1;
        p1_down[1] -= 1;
        temp_dist = l2distance( p1_down, p2 );
        dist = std::min( dist, temp_dist );

        Point p1_downright = p1;
        p1_downright[0] += 1;
        p1_downright[1] -= 1;
        temp_dist = l2distance( p1_downright, p2 );
        dist = std::min( dist, temp_dist );

        return dist;
    }

    void create_rips_filtration( const std::string& infile_points_name, PointContainer& points_values,
                                 IndexContainer& to_new_index, Pdsm::BoundaryMatrix& bdmat,
                                 Pdsm::Radii& radii_vector,
                                 std::vector< Rips::IndexContainer >& complex ) {
        Rips::read_points( infile_points_name, points_values );
        std::cout << "The number of points: " << points_values.size() << std::endl;

        Rips::DistanceContainer radii;
        std::vector< std::vector< std::pair< Rips::Index, Rips::Coeff>> > boundaries;

        const int64_t nr_vertices = points_values.size();
        const int64_t nr_edges = (nr_vertices * (nr_vertices - 1)) / 2;
        const int64_t nr_faces = (nr_vertices * (nr_vertices - 1) * (nr_vertices - 2)) / 6;
        boundaries.resize( nr_vertices + nr_edges + nr_faces );
        Rips::Index counter = -1;

        // 0 dim
        for ( int i = 0; i < nr_vertices; ++i ) {
            ++counter;
            Rips::IndexContainer temp = {i};
            complex.push_back( temp );
            radii.push_back( 0 );
        }

        std::map< Rips::IndexContainer, Rips::Index > vertices_to_index;
        // 1 dim
        for ( int i = 0; i < nr_vertices; ++i ) {
            for ( int j = i + 1; j < nr_vertices; ++j ) {
                ++counter;
                Rips::IndexContainer temp = {i, j};
                complex.push_back( temp );

                Rips::DistanceType radius = Rips::distance( points_values[i], points_values[j] ) / 2;
                radii.push_back( radius );

                boundaries[counter].push_back( {i, -1} );
                boundaries[counter].push_back( {j, 1} );
                vertices_to_index[{i, j}] = counter;
            }
        }


        // 2 dim
        for ( int i = 0; i < nr_vertices; ++i ) {
            for ( int j = i + 1; j < nr_vertices; ++j ) {
                for ( int k = j + 1; k < nr_vertices; ++k ) {
                    ++counter;
                    Rips::IndexContainer temp = {i, j, k};
                    complex.push_back( temp );

                    Rips::DistanceType radius1 = Rips::distance( points_values[i], points_values[j] ) / 2;
                    Rips::DistanceType radius2 = Rips::distance( points_values[i], points_values[k] ) / 2;
                    Rips::DistanceType radius3 = Rips::distance( points_values[j], points_values[k] ) / 2;
                    Rips::DistanceType radius = std::max( {radius1, radius2, radius3} );
                    radii.push_back( radius );

                    boundaries[counter].push_back( {vertices_to_index[{i, j}], 1} );
                    boundaries[counter].push_back( {vertices_to_index[{i, k}], -1} );
                    boundaries[counter].push_back( {vertices_to_index[{j, k}], 1} );
                }
            }
        }

        auto critical_radii = radii;
        std::sort( critical_radii.begin(), critical_radii.end());
        critical_radii.erase( std::unique( critical_radii.begin(), critical_radii.end()), critical_radii.end());

        Rips::IndexContainer index_map;
        std::vector< bool > visited( radii.size(), false );
        for ( auto critical_radius : critical_radii ) {
            for ( int idx = 0, size = radii.size(); idx < size; ++idx ) {
                if ( radii[idx] <= critical_radius && !visited[idx] ) {
                    index_map.push_back( idx );
                    visited[idx] = true;
                }
            }
        }

        to_new_index.resize( index_map.size());
        for ( Rips::Index idx = 0, size = index_map.size(); idx < size; ++idx ) {
            to_new_index[index_map[idx]] = idx;
        }

        // Check if point indices are invariant.
        for ( Rips::Index idx = 0; idx < nr_vertices; ++idx ) {
            if ( idx != to_new_index[idx] ) {
                std::cerr << "The function to_new_index ERROR!" << std::endl;
                throw;
            }
        }

        bdmat.resize_cols( boundaries.size());
        for ( Rips::Index idx = 0, size = boundaries.size(); idx < size; ++idx ) {
            Pdsm::Index new_idx = to_new_index[idx];
            if ( boundaries[idx].empty()) {
                bdmat.set_dim( new_idx, 0 );
            } else {
                bdmat.set_dim( new_idx, boundaries[idx].size() - 1 );
                for ( auto& index_coeff : boundaries[idx] ) {
                    bdmat.set_matrix_value( to_new_index[index_coeff.first], new_idx, Core::Zp( index_coeff.second ));
                }
            }
        }
        bdmat.prune();
        bool flag_print = (bdmat.cols() < 26);
        if ( flag_print ) {
            std::cout << "\nThe boundary matrix is:" << std::endl;
            bdmat.print_matrix();
        } else {
            std::cout << "The size of boundary matrix is: " << bdmat.cols() << std::endl;
        }

        Rips::DistanceContainer sorted_radii( radii.size());
        for ( Rips::Index idx = 0; idx < radii.size(); ++idx ) {
            sorted_radii[to_new_index[idx]] = radii[idx];
        }
        radii_vector = std::move( sorted_radii );
        std::cout << "nr of vertices: " << nr_vertices << std::endl;
        std::cout << "nr of edges: " << nr_edges << std::endl;
        std::cout << "nr of faces: " << nr_faces << std::endl;
    }

    void create_vertex_matrix( const IndexContainer& to_new_index, const std::vector< IndexContainer >& complex,
                               Pdsm::IndexMatrix& vertex_matrix_domain,
                               boost::unordered_map< Pdsm::IndexColumn, Pdsm::Index >& sorted_vertices_domain_to_index ) {
        vertex_matrix_domain.resize_cols( complex.size());
        for ( Rips::Index idx = 0, size = complex.size(); idx < size; ++idx ) {
            Pdsm::Index new_idx = to_new_index[idx];
            Pdsm::IndexColumn temp_vertices;
            for ( auto vertex : complex[idx] ) {
                temp_vertices.push_back( to_new_index[vertex] );
            }
            vertex_matrix_domain.set_col( new_idx, temp_vertices );
            auto sorted_vertices = temp_vertices;
            std::sort( sorted_vertices.begin(), sorted_vertices.end());
            sorted_vertices_domain_to_index[temp_vertices] = new_idx;
        }
    }

}

#endif //PDSM_RIPS_TORUS_H
