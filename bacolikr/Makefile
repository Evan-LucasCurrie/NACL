CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -Wpedantic -ffp-contract=off
AR ?= ar
SRC = src/abd.c src/linear.c src/spatial.c src/integrator.c src/solver.c
OBJ = $(patsubst src/%.c,build/obj/%.o,$(SRC))
.PHONY: all test clean
all: build/libbacolikr.a
build/libbacolikr.a: $(OBJ)
	$(AR) rcs $@ $^
build/obj/%.o: src/%.c src/abd.h src/linear.h src/spatial.h src/integrator.h src/bacolikr.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@
build/obj/spatial.o: src/coefficients.inc
build/tests/test_solver: tests/test_solver.c build/libbacolikr.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbacolikr.a -lm -o $@
build/tests/test_integrator: tests/test_integrator.c tests/test_data/integrator_fixtures.h tests/test_data/fixture_stream.h build/libbacolikr.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbacolikr.a -lm -o $@
build/tests/test_abd: tests/test_abd.c tests/test_data/abd_fixture.h tests/test_data/fixture_stream.h build/libbacolikr.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbacolikr.a -lm -o $@
build/tests/test_spatial: tests/test_spatial.c src/spatial.c src/coefficients.inc tests/test_data/spatial_fixture.h tests/test_data/fixture_stream.h build/libbacolikr.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbacolikr.a -lm -o $@
build/tests/test_features: tests/test_features.c build/libbacolikr.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbacolikr.a -lm -o $@
build/tests/test_integrator_optimization: tests/test_integrator_optimization.c build/libbacolikr.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbacolikr.a -lm -o $@
build/tests/test_linear_optimization: tests/test_linear_optimization.c build/libbacolikr.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbacolikr.a -lm -o $@
build/tests/test_solver_optimization: tests/test_solver_optimization.c build/libbacolikr.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbacolikr.a -lm -pthread -o $@
test: build/tests/test_solver build/tests/test_integrator build/tests/test_abd build/tests/test_spatial build/tests/test_features build/tests/test_integrator_optimization build/tests/test_linear_optimization build/tests/test_solver_optimization
	cd build && ./tests/test_abd
	cd build && ./tests/test_integrator
	cd build && ./tests/test_solver
	cd build && ./tests/test_spatial
	cd build && ./tests/test_features
	cd build && ./tests/test_integrator_optimization
	cd build && ./tests/test_linear_optimization
	cd build && ./tests/test_solver_optimization
clean:
	rm -rf build
