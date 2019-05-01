# Find module for RooUnfold -- by Dario Berzano <dario.berzano@cern.ch>
#
# From your `CMakeLists.txt` use it as simply as:
#
#   find_package(RooUnfold)
#
# It expects `-DROOUNFOLD=<prefix>` to be passed as argument, and it falls back
# to the environment variable `ROOUNFOLD_ROOT`.
#
# A target is exported if RooUnfold is found. The target takes care of defining
# all the library paths, compile flags, and include paths. Usage is very simple:
#
#   if(RooUnfold_FOUND)
#     target_link_library(MyAnalysisClass MyDep1 MyDep2 RooUnfold::RooUnfold)
#   endif()
#
# No additional command is required (no `include_directories()`, etc.) and no
# variable other than `RooUnfold_FOUND` is exported.

include(FindPackageHandleStandardArgs)

if(NOT DEFINED ROOUNFOLD)
  set(ROOUNFOLD "$ENV{ROOUNFOLD_ROOT}")
endif()

find_library(ROOUNFOLD_LIBPATH RooUnfold
             PATHS "${ROOUNFOLD}/lib" NO_DEFAULT_PATH)
find_path(ROOUNFOLD_INCLUDE_DIR RooUnfold.h
          PATHS "${ROOUNFOLD}/include" NO_DEFAULT_PATH)

find_package_handle_standard_args(RooUnfold DEFAULT_MSG
                                  ROOUNFOLD_LIBPATH ROOUNFOLD_INCLUDE_DIR)

# Set RooUnfold::RooUnfold target
add_library(RooUnfold::RooUnfold SHARED IMPORTED)
set_target_properties(RooUnfold::RooUnfold PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${ROOUNFOLD_INCLUDE_DIR}"
  IMPORTED_LOCATION             "${ROOUNFOLD_LIBPATH}"
  INTERFACE_COMPILE_DEFINITIONS "WITH_ROOUNFOLD")

# Unset RooUnfold variables
unset(ROOUNFOLD_LIBPATH)
unset(ROOUNFOLD_INCLUDE_DIR)
