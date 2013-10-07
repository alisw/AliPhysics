include("${Vc_CMAKE_MODULES_DIR}/AddTargetProperty.cmake")

include_directories("${Vc_INCLUDE_DIR}")

string(REPLACE "-Weffc++" "" CXXFLAGS "${CXXFLAGS}")
set(libvc_compile_flags "-DVC_COMPILE_LIB ${CXXFLAGS} -O3 ${Vc_DEFINITIONS}")
vc_compile_for_all_implementations(_objs src/trigonometric.cpp FLAGS ${libvc_compile_flags}
      ONLY SSE2 SSE3 SSSE3 SSE4_1 AVX SSE+XOP+FMA4 AVX+XOP+FMA4)
set(_srcs src/const.cpp src/cpuid.cpp src/support.cpp ${_objs})
if(USE_AVX)
   list(APPEND _srcs src/avx_sorthelper.cpp)
else()
   if(NOT Vc_AVX_INTRINSICS_BROKEN)
      # we'd still like to have avx_sorthelper.cpp built in, but that requires compilation with -mavx (or a comparable flag)
      foreach(_flag "-xAVX" "-mavx" "/arch:AVX")
         AddCompilerFlag("${_flag}" CXX_RESULT _flag_works)
         if(_flag_works)
            if(_flag STREQUAL "-xAVX")
               set(_flag "${_flag} -diag-disable 10121") # disable the warning "overriding -xSSE4.2 with -xAVX"
            endif()
            list(APPEND _srcs src/avx_sorthelper.cpp)
            set_source_files_properties(src/avx_sorthelper.cpp PROPERTIES COMPILE_FLAGS "${_flag}")
            break()
         endif()
      endforeach()
   endif()
endif()
add_library(Vc STATIC ${_srcs})
add_target_property(Vc COMPILE_FLAGS ${libvc_compile_flags})

# vim: ft=cmake sw=3 et
