# Source attribution

The algorithms and problem definitions are adapted from
[BACOLIKR (2024)](https://cs.smu.ca/~muir/BACOLI-3_Webpage.htm), distributed
through Saint Mary's University. They solve one-dimensional time-dependent
parabolic partial differential equations with event detection.

- Original BACOL: Copyright (c) 2002, Rong Wang, Pat Keast, Paul Muir.
- BACOLI: Copyright (c) 2019, Paul Muir, Jack Pew, Connor Tannahill.
- BACOLIKR wrapper and drivers: Copyright (c) 2024, Paul Muir, Jack Pew,
  Connor Tannahill.

The original BSD-style two-clause copyright notice and disclaimer are retained
in `LICENSE`. Time integration formulas follow the modified DASKR routines
distributed with BACOLIKR. Spline evaluation uses the Cox–de Boor recurrence;
linear systems use native C band and almost-block-diagonal elimination.
