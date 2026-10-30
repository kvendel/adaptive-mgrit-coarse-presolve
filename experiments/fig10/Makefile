XBRAID_DIR = ../../xbraid
BRAID_DIR = $(XBRAID_DIR)/braid
include $(XBRAID_DIR)/makefile.inc

BRAID_FLAGS = -I$(BRAID_DIR)
BRAID_LIB_FILE = $(BRAID_DIR)/libbraid.a

SRC_DIR = src

CC ?= gcc

SOURCES = $(wildcard $(SRC_DIR)/*.c)
SEQUENTIAL_SOURCES = $(foreach file,$(SOURCES),$(if $(findstring sequential,$(file)),$(file)))
PARALLEL_SOURCES = $(filter-out $(SEQUENTIAL_SOURCES),$(SOURCES))

SEQUENTIAL_TARGETS = $(SEQUENTIAL_SOURCES:.c=)
PARALLEL_TARGETS = $(PARALLEL_SOURCES:.c=)

TARGETS = $(SEQUENTIAL_TARGETS) $(PARALLEL_TARGETS)

.PHONY: all parallel sequential clean list xbraid

all: $(TARGETS)

parallel: $(PARALLEL_TARGETS)

sequential: $(SEQUENTIAL_TARGETS)

xbraid:
	$(MAKE) -C $(XBRAID_DIR)

$(PARALLEL_TARGETS): %: %.c $(BRAID_LIB_FILE)
	@echo "Building parallel executable $@ ..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $< -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

$(SEQUENTIAL_TARGETS): %: %.c
	@echo "Building sequential executable $@ ..."
	$(CC) $(CFLAGS) $< -o $@ $(LFLAGS)

list:
	@echo "Parallel sources:"
	@printf "  %s\n" $(PARALLEL_SOURCES)
	@echo "Sequential sources:"
	@printf "  %s\n" $(SEQUENTIAL_SOURCES)
	@echo "Parallel targets:"
	@printf "  %s\n" $(PARALLEL_TARGETS)
	@echo "Sequential targets:"
	@printf "  %s\n" $(SEQUENTIAL_TARGETS)

clean:
	rm -f $(TARGETS)