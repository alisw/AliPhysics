################################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
################################################################################

# Set the platform variable if not yet set.
ifndef PLATFORM
PLATFORM = $(shell uname)
endif

# Get verison numbers.
MAJOR_VERSION = $(shell cat $(SRC_DIR)/Version/major.number)
MINOR_VERSION = $(shell cat $(SRC_DIR)/Version/minor.number)


.SUFFIXES: .cxx .hpp .h .o .d

# Set the module name, i.e. the name of the executable / library.
ifdef BINARY
MODULE_NAME = $(BINARY)
endif

ifdef LIBRARY
MODULE_NAME = $(LIBRARY)
endif

################################################################################
# Directory layout:

TOP_DIR        = $(shell pwd)
BIN_DIR        = $(TOP_DIR)/bin
LIB_DIR        = $(TOP_DIR)/lib
INCLUDE_DIR    = $(TOP_DIR)/include
OUTPUT_DIR     = $(TOP_DIR)/output
BUILD_DIR      = $(TOP_DIR)/build
MACRO_DIR      = $(TOP_DIR)/macros
SRC_DIR        = $(TOP_DIR)/src

# The base target directories.
BIN_TARGET_DIR    = $(BIN_DIR)/$(PLATFORM)
LIB_TARGET_DIR    = $(LIB_DIR)/$(PLATFORM)
OUTPUT_TARGET_DIR = $(OUTPUT_DIR)/$(MODULE_NAME)/$(PLATFORM)

# The list of source code sub directory names.
SRC_SUBDIR_NAMES = AliRoot Buffers Clustering Control Decision PubSub PipeIO \
	Tracking Version System System/Linux System/Win32 Debug Framework DDL \
	BCMP

# Create full source and target directory names.
SRC_SUBDIRS    = $(addprefix $(SRC_DIR)/,$(SRC_SUBDIR_NAMES))
OUTPUT_TARGET_SUBDIRS = $(addprefix $(OUTPUT_TARGET_DIR)/,$(SRC_SUBDIR_NAMES))

# Setup the distribution file name and directory structure.
DISTRIBUTION_DIRNAME = dHLT_v$(MAJOR_VERSION)-$(MINOR_VERSION)
DISTRIBUTION_DIR = $(OUTPUT_DIR)/distrib/$(PLATFORM)
DISTRIBUTION_SUBDIR = $(DISTRIBUTION_DIR)/$(DISTRIBUTION_DIRNAME)

################################################################################

# Set the target name, symbolic link and full distribution directory paths.
ifdef BINARY
TARGET_NAME = $(BIN_TARGET_DIR)/$(BINARY)$(BIN_EXT)
SYMLINK_NAME = $(BIN_DIR)/$(BINARY)$(BIN_EXT)
DISTRIBUTION_TARGET_DIR = $(DISTRIBUTION_SUBDIR)/bin
endif

ifdef LIBRARY
TARGET_NAME = $(LIB_TARGET_DIR)/$(LIB_PREFIX)$(LIBRARY)$(LIB_EXT)
SYMLINK_NAME = $(LIB_DIR)/$(LIB_PREFIX)$(LIBRARY)$(LIB_EXT)
DISTRIBUTION_TARGET_DIR = $(DISTRIBUTION_SUBDIR)/lib
endif

################################################################################

include $(BUILD_DIR)/$(PLATFORM).platform.mk

################################################################################

# If we want a verbose build then all the actuall commands that are executed
# are printed to screen.
ifdef VERBOSE_MAKE
VERBOSE =
else
VERBOSE = @
endif

################################################################################
# Setup the CXXFLAGS compiler flags and LINKFLAGS for the linker.

ifdef BINARY
CXXFLAGS = $(INCLUDES) $(MACROS) $(BINARY_CXXFLAGS)
LINKFLAGS = $(BINARY_CXXFLAGS) $(LIBRARY_PATHS) $(LIBRARIES)
endif

ifdef LIBRARY
CXXFLAGS = $(INCLUDES) $(MACROS) $(LIBRARY_CXXFLAGS)
LINKFLAGS = $(LIBRARY_CXXFLAGS) $(LIBRARY_PATHS) $(LIBRARIES)
endif

ifdef DEBUGGING
MACROS += DEBUG
CXXFLAGS += $(CXX_DEBUG)
LINKFLAGS += $(CXX_DEBUG)
BIN_TARGET_DIR := $(BIN_TARGET_DIR)-debug
LIB_TARGET_DIR := $(LIB_TARGET_DIR)-debug
OUTPUT_TARGET_DIR := $(OUTPUT_TARGET_DIR)-debug
DISTRIBUTION_FILE = dHLT_$(PLATFORM)-debug_v$(MAJOR_VERSION)-$(MINOR_VERSION).tar.gz
else
CXXFLAGS += $(OPTIMIZATION)
LINKFLAGS += $(OPTIMIZATION)
DISTRIBUTION_FILE = dHLT_$(PLATFORM)_v$(MAJOR_VERSION)-$(MINOR_VERSION).tar.gz
endif

################################################################################

# Create a list of headers from the source file list.
HEADERS := $(SOURCES:.cxx=.hpp)

ifndef DICTIONARY_HEADERS
DICTIONARY_HEADERS = $(HEADERS)
endif

# Add the dictionary file without adding its header to the HEADERS list.
ifdef DICTIONARY
SOURCES += $(DICTIONARY)
endif

# Create a list of object .o and dependancy .d files from the list of source files.
OBJECTS = $(addprefix $(OUTPUT_TARGET_DIR)/,$(SOURCES:.cxx=.o))
DEPENDS = $(addprefix $(OUTPUT_TARGET_DIR)/,$(SOURCES:.cxx=.d))

# Make the includes search paths from the source directory and the output
# target directory.
INCLUDES += $(SRC_DIR) $(OUTPUT_TARGET_DIR)

# The library path should also contain the lib sub directory.
LIBRARY_PATHS += $(LIB_TARGET_DIR)

# Add compiler specific flag prefices:
INCLUDES := $(addprefix $(INCLUDE_PREFIX),$(INCLUDES))
MACROS := $(addprefix $(MACRO_PREFIX),$(MACROS))
LIBRARIES := $(addprefix $(LIBRARY_PREFIX),$(LIBRARIES))
LIBRARY_PATHS := $(addprefix $(LIBPATH_PREFIX),$(LIBRARY_PATHS))

################################################################################
# Tell make where to find source files.
# Note: this must come after the OUTPUT_TARGET_DIR is appended with -debug
vpath %.cxx $(SRC_DIR) $(SRC_SUBDIRS) $(OUTPUT_TARGET_DIR) $(OUTPUT_TARGET_SUBDIRS)
vpath %.hpp $(SRC_DIR) $(SRC_SUBDIRS) $(OUTPUT_TARGET_DIR) $(OUTPUT_TARGET_SUBDIRS)
vpath %.h $(SRC_DIR) $(SRC_SUBDIRS) $(OUTPUT_TARGET_DIR) $(OUTPUT_TARGET_SUBDIRS)
vpath %.number $(SRC_DIR)/Version $(OUTPUT_TARGET_DIR)
