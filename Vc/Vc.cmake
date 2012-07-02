set(Vc_INCLUDE_DIR          "${PROJECT_SOURCE_DIR}/Vc/include")
set(Vc_CMAKE_MODULES_DIR    "${PROJECT_SOURCE_DIR}/Vc/cmake")
include("${Vc_CMAKE_MODULES_DIR}/VcMacros.cmake")
vc_set_preferred_compiler_flags()

macro(ALICE_UseVc)
   include_directories(SYSTEM "${Vc_INCLUDE_DIR}")
   add_definitions(${Vc_DEFINITIONS})
endmacro()

# vim: ft=cmake sw=3 et
