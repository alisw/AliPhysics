include("${Vc_CMAKE_MODULES_DIR}/AddTargetProperty.cmake")

include_directories("${Vc_INCLUDE_DIR}")

string(REPLACE "-Weffc++" "" CXXFLAGS "${CXXFLAGS}")

set(libvc_compile_flags "-DVC_COMPILE_LIB ${CXXFLAGS} -O3")
vc_compile_for_all_implementations(_objs0 src/trigonometric.cpp FLAGS ${libvc_compile_flags}
      ONLY SSE2 SSE3 SSSE3 SSE4_1 AVX SSE+XOP+FMA4 AVX+XOP+FMA4)
vc_compile_for_all_implementations(_objs1 src/avx_sorthelper.cpp FLAGS ${libvc_compile_flags} ONLY AVX)
set(_srcs src/const.cpp src/cpuid.cpp src/support.cpp ${_objs0} ${_objs1})
add_library(Vc STATIC ${_srcs})
add_target_property(Vc COMPILE_FLAGS "${libvc_compile_flags} ${Vc_DEFINITIONS}")

# vim: ft=cmake sw=3 et
