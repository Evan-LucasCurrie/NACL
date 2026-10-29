# COLNEWSC in C

Standalone adaptive Gaussian collocation for first-order ordinary differential equation boundary value problems. The library provides a primary solver and a separate reference solver with upstream-compatible iteration and interpolation policies.

```sh
make             # build/libcolnewsc.a
make test        # numerical, API, ownership and exact-reference tests
```

Implementation and headers live in `src/`, and tests in `tests/`.
Generated objects, libraries, executables and debug files stay
in `build/`; `make clean` removes that directory. Compile callers with `-Isrc`.

Include `colnewsc.h` for the primary solver or `colnewsc_reference.h` for the reference solver, and link `build/libcolnewsc.a -lm`.
Callback matrices are row major. The primary solver
supports 1-160 first-order equations; the reference solver supports 1-20.

Create options with `colnewsc_default_options()`, supply callbacks and interval bounds in `colnewsc_problem`, call `colnewsc_solve`, sample the continuous solution with `colnewsc_evaluate`, and release it with `colnewsc_free`. An optional initial solution supports parameter continuation. Returned primary solutions own their numerical state and can outlive their callbacks.

Reference solutions retain the problem and callback context and re-evaluate callbacks during interpolation; those objects must outlive the solution. Separate solutions have independent state. Shared evaluation of one reference solution requires external locking.

## Numerical methods

The primary solver includes Gauss orders 1-7, integrated Lagrange bases,
condensed-stage systems, affine-invariant damped Newton iteration, analytic or centered numerical derivatives, Richardson error estimates, adaptive meshes, fixed interior mesh points, and continuous superconvergent output for orders 2-4. Other orders use ordinary integrated-basis output.

The reference solver retains one-based, column-major workspaces and supported aliases, original mesh selection and doubling, frozen Jacobians, block elimination, LINPACK normal/transposed solves, streaming BLAS operations, Vandermonde construction, and all continuous interpolation branches.

Stage and global LU storage is reused between Newton iterations. Small stage systems use shared factor slabs; larger uncoupled systems retain compact band storage. Specialized scalar and two-component kernels retain accumulation order, zero-weight skips and alias fallbacks. Floating-point contraction is disabled and no fast-math flags are used.

## Tests

Test sources, analytic equations and frozen fixtures live in `tests/`.

Tests cover every Gauss order, analytic and numerical derivatives, nonlinear Bratu and sine equations, continuous polynomial solutions, custom and fixed meshes, continuation, pivots, dense/band transitions, buffer reuse, ownership, invalid callbacks and interpolation diagnostics. The reference suite checks **704 sampled values, 24,275 stored coefficients, and every mesh history** against frozen records embedded as C data in `tests/test_data/reference_fixtures.h`. It also tests 84 order/dimension/boundary combinations. No external fixture files or sibling projects are needed.

## Attribution

COLNEWSC is credited to Mark Adams and Paul Muir as an extension of COLNEW. The original numerical coefficients and algorithm operation order are retained.
