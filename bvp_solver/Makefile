CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -pedantic -ffp-contract=off
LDLIBS ?= -lm
OBJECTS = build/obj/method.o build/obj/linalg.o build/obj/solver.o
.PHONY: all test clean
all: build/libbvp_solver.a
build/libbvp_solver.a: $(OBJECTS)
	ar rcs $@ $^
build/obj/%.o: src/%.c src/bvp.h src/internal.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@
build/tests/tests: tests/tests.c tests/test_problems.h build/libbvp_solver.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbvp_solver.a $(LDLIBS) -o $@
build/tests/test_optimization: tests/test_optimization.c build/libbvp_solver.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libbvp_solver.a $(LDLIBS) -o $@
test: build/tests/tests build/tests/test_optimization
	cd build && ./tests/tests
	cd build && ./tests/test_optimization
clean:
	rm -rf build
