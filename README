# pdsm (persistence diagrams of sampled maps)

This program is an implementation for the paper "The persistent homology of a sampled map: From a viewpoint of quiver representations" https://arxiv.org/abs/1810.11774. This program calculates the birth-death pairs in the persistence diagram of a sampled map using Vietoris–Rips complexes.

## Requirements

Boost https://www.boost.org/

Eigen http://eigen.tuxfamily.org/

## Remarks on the source codes

The arXiv peper is v2 but the source code still use the terminology for v1. The author will update it soon.

## Build

```shell
mkdir build
cd build
cmake ..
make
```

or

```shell
mkdir build
cd build
g++ ../main.cpp -o pdsm -std=c++1z -lstdc++fs -I/usr/local/include/eigen3 -O3
g++ ../main_torus.cpp -o pdsm_torus -std=c++1z -lstdc++fs -I/usr/local/include/eigen3 -O3
```

## Usage

The directory ./data stores some sample data sets for pdsm.

### Example 1 (twice mapping on a circle w/ 20 sample points w/o noise)

```shell
./pdsm ../data/points/circle_points20_2_0.txt ../data/points/circle_points_codomain20_2_0.txt ../data/map/map20.txt ../data/output20_2_0
```

- ../data/points/circle_points20_2_0.txt - a point cloud for domain circle.

- ../data/points/circle_points_codomain20_2_0.txt - a point cloud for codomain circle.

- ../data/map/map20.txt - (trivial) mapping list (the i-th line has the number i, so this maps the i-th point in circle_points20_2_0.txt to the i-th point in circle_points_codomain20_2_0.txt).

- ../data/output20_2_0 - results are written in this directory.

  The directory output20_2_0 will store two .txt files: bd0.txt and generator0.txt.

  - bd0.txt - the unique birth-death pair values.
  - generator0.txt - the corresponding generator in domain, where number i shows the i-th point in circle_points20_2_0.txt.

### Example 2 (twice mapping on a circle w/ 100 sample points w/ Gaussian noise of \sigma = 0, 0.09, 0.18, 0.3)

```shell
./pdsm ../data/points/circle_points100_2_0.txt ../data/points/circle_points_codomain100_2_0.txt ../data/map/map100.txt ../data/output100_2_0
./pdsm ../data/points/circle_points100_2_0.09.txt ../data/points/circle_points_codomain100_2_0.09.txt ../data/map/map100.txt ../data/output100_2_0.09
./pdsm ../data/points/circle_points100_2_0.18.txt ../data/points/circle_points_codomain100_2_0.18.txt ../data/map/map100.txt ../data/output100_2_0.18
./pdsm ../data/points/circle_points100_2_0.3.txt ../data/points/circle_points_codomain100_2_0.3.txt ../data/map/map100.txt ../data/output100_2_0.3
```

### Example 3 (mapping on a torus w/ 8^2 sample points w/o noise)

```shell
./pdsm_torus ../data/points/torus_points64_2111_0.txt ../data/points/torus_points_codomain64_2111_0.txt ../data/map/map64.txt ../data/output_torus_64_2111_0
```

## Acknowledgements

gyoza https://bitbucket.org/remere/gyoza developed by Emerson G. Escolar

PHAT https://github.com/blazs/phat developed by Ulrich Bauer, Michael Kerber, Jan Reininghaus

