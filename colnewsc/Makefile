CC ?= cc
CFLAGS ?= -std=c11 -O3 -Wall -Wextra -Wpedantic -ffp-contract=off
AR ?= ar
LDLIBS ?= -lm
SOURCES = src/colnewsc.c src/basis.c src/interpolant.c src/linalg.c src/reference/algebra.c src/reference/interpolation.c src/reference/mesh.c src/reference/control.c src/reference/support.c src/reference/context.c
OBJECTS = $(patsubst src/%.c,build/obj/%.o,$(SOURCES))
.PHONY: all test clean
all: build/libcolnewsc.a
build/libcolnewsc.a: $(OBJECTS)
	$(AR) rcs $@ $^
build/obj/%.o: src/%.c src/colnewsc.h src/internal.h src/colnewsc_reference.h src/reference/reference_internal.h src/reference/methods.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@
build/tests/tests: tests/tests.c tests/problems.c tests/problems.h tests/workloads.h build/libcolnewsc.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< tests/problems.c build/libcolnewsc.a $(LDLIBS) -o $@
build/tests/reference_tests: tests/reference_tests.c tests/problems.c tests/problems.h tests/reference_problems.c tests/test_data/reference_fixtures.h build/libcolnewsc.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< tests/problems.c tests/reference_problems.c build/libcolnewsc.a $(LDLIBS) -o $@
build/tests/primary_optimization_tests: tests/primary_optimization_tests.c build/libcolnewsc.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libcolnewsc.a $(LDLIBS) -o $@
build/tests/reference_optimization_tests: tests/reference_optimization_tests.c build/libcolnewsc.a
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -Isrc $< build/libcolnewsc.a $(LDLIBS) -o $@
test: build/tests/tests build/tests/reference_tests build/tests/primary_optimization_tests build/tests/reference_optimization_tests
	cd build && ./tests/tests
	cd build && ./tests/reference_tests
	cd build && ./tests/primary_optimization_tests
	cd build && ./tests/reference_optimization_tests
clean:
	rm -rf build
