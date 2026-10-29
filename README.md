# Numerical analysis

Three standalone numerical libraries, with C tests.

| Library    | Methods                                                                                                                                                          | Public headers                                                                   |
| ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| COLNEWSC   | Adaptive Gaussian collocation, orders 1–7; Newton iteration, continuous superconvergent interpolation and continuation; an additional upstream-compatible solver | `colnewsc/src/colnewsc.h`, `colnewsc/src/colnewsc_reference.h`                   |
| BVP_SOLVER | MIRK orders 2/4/6, adaptive meshes, error estimators, unknown parameters, singular endpoints and continuous output                                               | `bvp_solver/src/bvp.h`                                                           |
| BACOLIKR   | C1 spline collocation, orders 2–10, LOI/SCI spatial estimates, BDF orders 1–5, adaptive remeshing, conservative fluxes, DAE initialization and events            | `bacolikr/src/bacolikr.h`, `bacolikr/src/integrator.h`, `bacolikr/src/spatial.h` |

These were all ported from Fortran to C for portability purposes.

Requirements: a C compiler, the Make tool, and the system math library. The concurrent cache test uses POSIX threads. Each library has its own Makefile and can be built independently.

Each library uses the same directory layout:

```text
src/          solver implementation, and headers
tests/        test sources, fixtures and test executables
Makefile      builds the library, tests
```

Static libraries are built in each library's root directory. Use that library's
`src/` directory as the include path when compiling a caller.

From this directory:

```sh
make             # build the three static libraries
make test        # run numerical and API regression tests
make sanitize    # run C tests with address/undefined-behavior sanitizers
make clean       # remove generated objects, libraries and executables
```

Builds use double precision, optimization level 3 and floating-point contraction disabled. No fast-math flags are used. Override `CC` and `CFLAGS` as needed.

See the library READMEs for callback layouts and ownership. Original copyright
and third-party notices are retained in the library folders.
