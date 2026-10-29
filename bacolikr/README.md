# BACOLIKR in C

A standalone C library for systems of one-dimensional parabolic PDEs. It uses
spline collocation, adaptive BDF integration, spatial error control and event
detection.

```sh
make
make test
```

Implementation and headers live in `src/`, and tests in `tests/`.
Generated objects, libraries, executables and debug files stay
in `build/` `make clean` removes that directory. Compile callers with `-Isrc`,
include `bacolikr.h`, and link `build/libbacolikr.a` with `-lm`. Supply a `BacProblem`
with component count, context and callbacks, call `bac_solver_new`, then
`bac_solver_advance`. `BacView` evaluates the current spline and performs Gaussian quadrature. Events receive a view of the dense solution. Release successful or partially initialized solvers with `bac_solver_free`, and advance results with bac_advance_result_free`. The caller owns callback contexts. Endpoint `side` is zero for left and one for right. Jacobian callbacks return one when supplied and zero to request numerical differences. Arrays are component interleaved, derivative matrices are row major.

The implementation includes spline orders 2-10, LOI/SCI spatial estimates,
BDF orders 1-5, verified and DASKR time controllers, colored differences,
consistent DAE initialization, dense output, directional roots, conservative
fluxes and warm/cold remeshing. Solver spatial adaptation uses orders 3-10.

Optimizations share immutable Newton factors across snapshots, reuse private
workspaces and snapshot vectors, and retain active LU extents. Accepted grids
are borrowed until remeshing requires ownership transfer. Quadrature and view
cache buffers are built on demand. Scalar spline sums and conservative assembly
avoid redundant work. Numerical methods, tolerances and stopping criteria are
unchanged. Floating-point contraction is disabled and fast-math is not used.

Tests include analytic and nonlinear solutions, all spline orders, both spatial
estimators, multi-component systems, Neumann/Robin boundaries, conservative
Jacobians, mass conservation, event restart, failed-advance rollback and cache
ownership. Test sources and helpers live in `tests/`. Frozen reference records
are embedded as C headers in `tests/test_data/`

Copyright and upstream attribution are retained in `LICENSE` and
`THIRD_PARTY_NOTICES.md`.
