# MMA and GCMMA
*A self-contained C++ implementation of MMA and GCMMA.*

[![Build Status](https://travis-ci.org/jdumas/mma.svg?branch=master)](https://travis-ci.org/jdumas/mma)

This repository contains single-file C++ implementations of MMA and GCMMA, as described in [1,2,3].
The code in this repository is based on an original code by Niels Aage, using the subproblem solver described in [4].
I made some modifications to the original code to be more C++-like (e.g. using `std::vector<>`), and slightly changed the way bounds are computed in the subproblem approximation.
I have also extended the solver to use the GCMMA method [2,3]. Both the MMA and GCMMA solvers can be found in their respective folder in `src/`.

The goal of this project is to provide a minimal, self-contained implementation of MMA and GCMMA which are easy to use and easy to integrate in an existing C++ project.
The code is still in **beta**, meaning it hasn't been tested extensively, so use at your own risks.

## Usage

This project use CMake, but you should be able to copy the files directly into your own project, since there are no dependencies.
Usage should be fairly obvious looking at the header files, but you can also take a look at the examples under the `tests/` folder.

## Other Projects

- [TopOpt_in_PETSc](https://github.com/topopt/TopOpt_in_PETSc): Multi-processor implementation of MMA using MPI and PETCs.
- [NLopt](https://nlopt.readthedocs.io/en/latest/): Open-source library for nonlinear optimization with bindings in different languages. Contains a different implementation of GCMMA.

## References

1. Svanberg, K. (1987). “The Method of Moving Asymptotes—a New Method for Structural Optimization”. International Journal for Numerical Methods in Engineering, 24(2), 359–373. https://doi.org/10.1002/nme.1620240207.
2. Svanberg, K. (2002). “A Class of Globally Convergent Optimization Methods Based on Conservative Convex Separable Approximations”. SIAM Journal on Optimization, 12(2), 555–573. https://doi.org/10.1137/s1052623499362822.
3. Svanberg, K. (2007). MMA and GCMMA - Two Methods for Nonlinear Optimization. Technical report. https://people.kth.se/~krille/mmagcmma.pdf
4. Aage, N., & Lazarov, B. S. (2013). Parallel framework for topology optimization using the method of moving asymptotes. Structural and Multidisciplinary Optimization, 47(4), 493–505. https://doi.org/10.1007/s00158-012-0869-2
