# SPDX-License-Identifier: MIT-0
# Standalone GNU Make build; does not include or invoke the parent makefile.
# Expected submodule location: VERKKO_ROOT/src/hpc-consensus.
.DEFAULT_GOAL := all
VERKKO_ROOT ?= $(abspath ../..)
VERKKO_SRC ?= $(VERKKO_ROOT)/src
VERKKO_BUILD ?= $(VERKKO_ROOT)/build
VERKKO_LIB ?= $(VERKKO_BUILD)/lib/libverkko.a
BUILD_DIR ?= build
PKG_CONFIG ?= pkg-config
PYTHON ?= python3

# Use Make's standard CXX default or the user's environment/command-line setting.
# Compiler selection belongs to the caller, including the parent Verkko build.
ifneq ($(shell uname -s),Darwin)
PLATFORM_LIBS := -ldl
endif
CXXFLAGS ?= -O3 -g
# GCC 9 accepts the draft spelling, not -std=c++20. Newer GCC and Clang
# accept c++2a too; leave compiler selection entirely with the caller.
CXX_STANDARD ?= c++2a
OPENMP_FLAGS ?= -fopenmp
OPENMP_LIBS ?= $(OPENMP_FLAGS)
DEP_CFLAGS := $(shell $(PKG_CONFIG) --cflags openssl libcurl liblzma 2>/dev/null)
DEP_LIBS := $(shell $(PKG_CONFIG) --libs openssl libcurl liblzma 2>/dev/null)
# Quote-only search avoids utility/math.H shadowing the system math.h on macOS.
PROJECT_FLAGS = -std=$(CXX_STANDARD) -pthread $(OPENMP_FLAGS) -Isrc \
  -iquote$(VERKKO_SRC)/canu/src/utility/src -iquote$(VERKKO_SRC)/canu/src/stores \
  $(DEP_CFLAGS)
PROJECT_LIBS = $(VERKKO_LIB) $(DEP_LIBS) -lz -lbz2 -pthread $(OPENMP_LIBS) $(PLATFORM_LIBS)
SOURCES := core models polish ends offsets hp_refinement io options
OBJECTS := $(addprefix $(BUILD_DIR)/obj/,$(addsuffix .o,$(SOURCES)))
PROGRAM := $(BUILD_DIR)/bin/hpc-consensus
TEST_NAMES := core offsets maintenance hp_refinement hp_flanks linked_events
TEST_BINS := $(addprefix $(BUILD_DIR)/tests/,$(TEST_NAMES))
HELPERS := $(BUILD_DIR)/bin/extract-regions $(BUILD_DIR)/bin/apply-run-model

.PHONY: all check-deps test clean
all: $(PROGRAM)

# Order-only prerequisite validates configuration without forcing recompilation.
check-deps:
	@test -f "$(VERKKO_SRC)/canu/src/stores/tgTig.H" || { echo 'Set VERKKO_ROOT to a Verkko checkout (expected ../../src/canu).' >&2; exit 1; }
	@test -f "$(VERKKO_LIB)" || { echo 'Build Verkko first; missing $(VERKKO_LIB)' >&2; exit 1; }
	@$(PKG_CONFIG) --exists openssl libcurl liblzma || { echo 'pkg-config requires openssl, libcurl and liblzma development files.' >&2; exit 1; }

$(BUILD_DIR)/obj/%.o: src/%.cpp Makefile | check-deps
	@mkdir -p "$(@D)"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) -MMD -MP -c $< -o $@

$(PROGRAM): $(OBJECTS) $(BUILD_DIR)/obj/main.o $(VERKKO_LIB)
	@mkdir -p "$(@D)"
	$(CXX) $(LDFLAGS) $(OBJECTS) $(BUILD_DIR)/obj/main.o $(PROJECT_LIBS) $(LDLIBS) -o $@

$(BUILD_DIR)/tests/%: tests/%.cpp $(OBJECTS) $(VERKKO_LIB) Makefile | check-deps
	@mkdir -p "$(@D)"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) -MMD -MP -MF $@.d $< $(OBJECTS) $(LDFLAGS) $(PROJECT_LIBS) $(LDLIBS) -o $@

$(VERKKO_LIB):
	@echo 'Build Verkko first; missing $@. Override VERKKO_ROOT or VERKKO_BUILD if needed.' >&2
	@exit 1

$(BUILD_DIR)/bin/extract-regions: tests/extract_regions.cpp $(OBJECTS) $(VERKKO_LIB) Makefile | check-deps
	@mkdir -p "$(@D)"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $< $(OBJECTS) $(LDFLAGS) $(PROJECT_LIBS) $(LDLIBS) -o $@

$(BUILD_DIR)/bin/apply-run-model: tests/apply_run_model.cpp $(OBJECTS) $(VERKKO_LIB) Makefile | check-deps
	@mkdir -p "$(@D)"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $< $(OBJECTS) $(LDFLAGS) $(PROJECT_LIBS) $(LDLIBS) -o $@

test: $(PROGRAM) $(TEST_BINS) $(HELPERS)
	@set -e; for t in $(TEST_BINS); do echo "Testing $$t"; "$$t"; done
	$(PYTHON) tests/hp_refinement_native.py $(abspath $(PROGRAM)) $(VERKKO_ROOT)
	$(PYTHON) tests/end_flank_native.py $(abspath $(PROGRAM)) $(VERKKO_ROOT)
	$(PYTHON) tests/offset_native.py $(abspath $(PROGRAM)) $(VERKKO_ROOT)
	$(PYTHON) tests/integration.py $(abspath $(PROGRAM)) $(VERKKO_ROOT)
	$(PYTHON) tests/calibration.py

clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.d) $(BUILD_DIR)/obj/main.o $(BUILD_DIR)/obj/main.d $(PROGRAM) $(TEST_BINS) $(addsuffix .d,$(TEST_BINS)) $(HELPERS)

-include $(OBJECTS:.o=.d) $(BUILD_DIR)/obj/main.d $(addsuffix .d,$(TEST_BINS))
