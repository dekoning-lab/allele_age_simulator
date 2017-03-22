# Allele age simulator

This is a simulator of reverse markov process that can be used to simulate a distribution of allele ages.

The process of interes is a Wright-Fisher population model, where the starting allele frequency is `1`, and observed allele frequency is `x`.

The simulator first creates a reversed process of the given chain, and then simulates from it.

The reverse Markov process is implemented according to [Chae and Kim, 1994](http://www.sciencedirect.com/science/article/pii/0167637794900205)

## Dependencies

This project depends on `armadillo`, `gsl`, and `openmp`.

## Installation

With the dependencies installed:

```
mkdir build && cd build
cmake ..
make

./allele_age_simulator --help
```

### Mac OS notes

This code was tested with `GCC`, since it fully support openmp.
Therefore, it is necessary to export the paths to the `GCC` executables before running `cmake`.

```
mkdir build && cd build

brew install gcc openmp armadillo gsl

export CC=/usr/local/bin/gcc-6
export CXX=/usr/local/bin/g++-6

cmake ..
make
```

### Linux notes

This code was tested against [linuxbrew](http://linuxbrew.sh/) installation.
Therefore, with `linuxbrew` in place, the compilation process should be the same as Mac OS.
