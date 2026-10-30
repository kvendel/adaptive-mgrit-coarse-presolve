# Top-level build configuration
#
# Default compiler settings target the Komondor Cray programming environment.
# They can be overridden from the command line, for example:
#   make CC=gcc MPICC=mpicc MPICXX=mpicxx

FIGURES := 6 7 8 9 10 11 12 13 14-16
FIGURE_TARGETS := $(addprefix fig,$(FIGURES))

CC ?= cc
MPICC ?= cc
MPICXX ?= CC

XBRAID_DIR := xbraid
XBRAID_BRAID_DIR := $(XBRAID_DIR)/braid
XBRAID_LIB := $(XBRAID_BRAID_DIR)/libbraid.a

.PHONY: all xbraid experiments $(FIGURE_TARGETS) list clean clean-experiments clean-xbraid

all: experiments

# User-facing target for building only the XBraid library.
xbraid: $(XBRAID_LIB)

# Build only libbraid.a, not the XBraid examples.
$(XBRAID_LIB):
	@echo "Building the XBraid library with CC=$(CC), MPICC=$(MPICC), MPICXX=$(MPICXX)..."
	$(MAKE) -C $(XBRAID_BRAID_DIR) CC=$(CC) MPICC=$(MPICC) MPICXX=$(MPICXX)

experiments: $(FIGURE_TARGETS)

# Every experiment requires libbraid.a, but the library is built only when needed.
$(FIGURE_TARGETS): $(XBRAID_LIB)
	@echo "Building programs for $@..."
	$(MAKE) -C experiments/$@ CC=$(CC) MPICC=$(MPICC) MPICXX=$(MPICXX)

list:
	@for fig in $(FIGURE_TARGETS); do \
		echo ""; \
		echo "Sources for $$fig:"; \
		$(MAKE) -C experiments/$$fig list || exit $$?; \
	done

clean: clean-experiments clean-xbraid

clean-experiments:
	@for fig in $(FIGURE_TARGETS); do \
		echo "Cleaning $$fig..."; \
		$(MAKE) -C experiments/$$fig clean || exit $$?; \
	done

clean-xbraid:
	@echo "Cleaning the XBraid library..."
	$(MAKE) -C $(XBRAID_BRAID_DIR) clean
