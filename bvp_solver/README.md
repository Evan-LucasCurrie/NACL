# BVP_SOLVER in C

Standalone numerical library for ordinary differential equation boundary-value
problems.

```sh
make
make test
```

Implementation and headers live in `src/`, and tests in `tests/`.
Generated objects, libraries, executables and debug files stay
in `build/`; `make clean` removes that directory. Compile callers with `-Isrc`,
include `bvp.h`, and link `build/libbvp_solver.a -lm`. Supply row-major callbacks in a `bvp_problem`, create an initial guess, call `bvp_solve`, and sample the returned continuous solution with `bvp_evaluate`. Release allocations with
`bvp_solution_free` and `bvp_guess_free`. Optional Jacobian callbacks may be null to request forward differences. Callback contexts and initial guesses remain
owned by the caller; solutions own their arrays and retain no callbacks.
Independent solver instances share no mutable global state.

## Algorithms

- MIRK orders 2, 4 and 6, with continuous value and derivative polynomials and
  an order-eight error-estimation tableau.
- Adaptive mesh redistribution, coarsening and retry; five error-control modes,
  defect sampling, deferred correction, Richardson estimates and condition
  estimates.
- Unknown parameters, singular-endpoint projections, damped Newton iteration,
  contracting chord steps and analytic or numerical derivatives.
- COLROW factorization and normal/transpose substitution, generic pivoted band
  LU, and cyclic one-sided Jacobi SVD pseudoinverses.
- Continuous interpolation, continuation, interval extension and versioned
  solution persistence. Saved files start with `BVP_SOLVER_C 1`.

Create options with `bvp_default_options`. The public header defines method,
tolerance, mesh-limit, error-control and diagnostic options. Unavailable optional
diagnostics are `NAN`; failures return `bvp_status` values.

## Tests

Test sources and shared problem definitions live in `tests/`.

`make test` covers 45 solver combinations across all MIRK orders and error
controls, supplied and numerical Jacobians, unknown parameters, singular
Lane-Emden, nonlinear Cash21 and Bratu equations, normal/transpose linear solves, rank-deficient pseudoinverses, continuation, extension, persistence and invalid inputs. Additional checks exercise zeroed scratch reuse, detached solution
lifetimes, cache bounds, nonlocal error cleanup and 17-state heap workspaces.

Temporary arrays reuse a per-solve cache capped at 1 MiB of retained buffer
capacity. Common small systems use bounded stack workspaces, and mesh transfer walks sorted nodes in order. Reused storage is zeroed before use; successful and failed solves release all cached storage.

Copyright notices and redistribution terms are retained in `LICENSE` and
`THIRD_PARTY_NOTICES.md`.
