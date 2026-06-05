# Build settings for i-tiger
SRC_DIR ?= src
BUILD_DIR ?= ./build
BIN_DIR ?= ./bin
TARGET ?= $(BIN_DIR)/i-tiger.x
TMP_FILE ?= /tmp/makedepf90_$(USER)

MPIF90 ?= mpiifort
MPI_INCLUDE ?= $(I_MPI_ROOT)/include/
LD := $(MPIF90)
LDFLAGS ?= 
LDLIBS ?= $(BLAS_LIBS)

BLAS_FLAGS ?= /opt/intel/oneapi/mkl/2022.1.0/include/
BLAS_LIBS ?= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

MOD_FLAG ?= -I
MODFLAGS ?= $(MOD_FLAG)$(BUILD_DIR) \
  $(MOD_FLAG)$(BLAS_FLAGS) \
  $(MOD_FLAG)$(MPI_INCLUDE)
D__FLAGS := -D__MPI -D__DFTI


FFLAGS_COMMON = -assume byterecl -g -traceback \
	-no-wrap-margin -nomodule -fpp -allow nofpp_comments \
  -stand f18 \
	$(MODFLAGS) $(D__FLAGS) -module $(BUILD_DIR)
FFLAGS_OPT ?= -O2 $(FFLAGS_COMMON)
FFLAGS_DEBUG ?= -O0 -check all -fpe0 $(FFLAGS_COMMON)

BUILD ?= release
ifeq ($(BUILD),debug)
  FFLAGS := $(FFLAGS_DEBUG)
else ifeq ($(BUILD),release)
  FFLAGS := $(FFLAGS_OPT)
else
  $(error Unsupported BUILD='$(BUILD)'. Use BUILD=debug or BUILD=release)
endif

# Collect all Fortran sources under src (recursive)
SRC := $(sort $(shell find $(SRC_DIR) -type f -name '*.f90' 2>/dev/null))
SRC_BASENAMES := $(notdir $(SRC))
SRC_DIRS := $(sort $(dir $(SRC)))

# This build layout assumes unique basenames for *.f90 files.
ifneq ($(words $(SRC_BASENAMES)),$(words $(sort $(SRC_BASENAMES))))
$(error Duplicate source basenames detected under $(SRC_DIR). Rename files to unique basenames.)
endif

# Objects are placed in build/ with basename mapping.
OBJ := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(notdir $(SRC))))
# VPATH lets pattern rules resolve basename-only prerequisites from src subdirs.
VPATH := $(SRC_DIRS)
DEP_FILE ?= depend.mk
