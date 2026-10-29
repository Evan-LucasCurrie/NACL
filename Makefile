LIBRARIES = colnewsc bvp_solver bacolikr
SANITIZER_CFLAGS = -std=c11 -O1 -g -Wall -Wextra -Wpedantic -ffp-contract=off -fsanitize=address,undefined -fno-omit-frame-pointer
.PHONY: all test sanitize clean
all test clean:
	@set -e; for library in $(LIBRARIES); do $(MAKE) -C $$library $@; done
sanitize:
	@set -e; for library in $(LIBRARIES); do \
		$(MAKE) -C $$library clean; \
		$(MAKE) -C $$library test CFLAGS="$(SANITIZER_CFLAGS)"; \
		$(MAKE) -C $$library clean; \
	done
