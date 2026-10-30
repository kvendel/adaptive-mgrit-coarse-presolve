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

.PHONY: all xbraid experiments $(FIGURE_TARGETS) list clean clean-experiments clean-xbraid

all: xbraid experiments

xbraid:
	@echo "Building XBraid with CC=$(CC), MPICC=$(MPICC), MPICXX=$(MPICXX)..."
	$(MAKE) -C xbraid CC=$(CC) MPICC=$(MPICC) MPICXX=$(MPICXX)

experiments: $(FIGURE_TARGETS)

$(FIGURE_TARGETS): xbraid
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
	@echo "Cleaning XBraid..."
	$(MAKE) -C xbraid clean
