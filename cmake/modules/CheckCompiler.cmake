# -*- mode: cmake -*-
MACRO ( Check_Compiler )


# Set a default build type for single-configuration
# CMake generators if no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

if (CMAKE_SYSTEM_NAME MATCHES Linux)
   MESSAGE("--- Found a Linux ssytem")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")
#      set ( CMAKE_SHARED_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_SHARED_LINKER_FLAGS}")
#      set ( CMAKE_MODULE_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_MODULE_LINKER_FLAGS}")
      # we profile...
#      if(CMAKE_BUILD_TYPE_TOLOWER MATCHES profile)
#        set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
#        set (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
#      endif(CMAKE_BUILD_TYPE_TOLOWER MATCHES profile)

   # Select flags.
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
   set(CMAKE_CXX_FLAGS_RELEASE        "-O2")
#   set(CMAKE_CXX_FLAGS_DEBUG          "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
   set(CMAKE_CXX_FLAGS_DEBUG          "-g")
   set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3 -fno-inline")
   set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
   set(CMAKE_CXX_FLAGS_ARRAY_CHECK    "-g3 -fno-inline -ftest-coverage -fprofile-arcs -fstack-protector")
   set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")
   set(CMAKE_C_FLAGS_RELEASE          "-O2")
#   set(CMAKE_C_FLAGS_DEBUG            "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
   set(CMAKE_C_FLAGS_DEBUG            "-g")
   set(CMAKE_C_FLAGS_DEBUGFULL        "-g3 -fno-inline")
   set(CMAKE_C_FLAGS_PROFILE          "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
   set(CMAKE_C_FLAGS_ARRAY_CHECK      "-g3 -fno-inline -ftest-coverage -fprofile-arcs -fstack-protector")

#   set ( CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
#   set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
   endif (CMAKE_COMPILER_IS_GNUCXX)

   if (CMAKE_C_COMPILER MATCHES "icc")
      MESSAGE("--- Found Intel compiler collection")
#      set ( CMAKE_SHARED_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_SHARED_LINKER_FLAGS}")
#      set ( CMAKE_MODULE_LINKER_FLAGS "-Wl,--fatal-warnings -Wl,--no-undefined -lc ${CMAKE_MODULE_LINKER_FLAGS}")

   # Select flags.
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
   set(CMAKE_CXX_FLAGS_RELEASE        "-O2")
   set(CMAKE_CXX_FLAGS_DEBUG          "-O2 -g -0b0 -noalign")
   set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g -Ob0 -noalign")
   set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")
   set(CMAKE_C_FLAGS_RELEASE          "-O2")
   set(CMAKE_C_FLAGS_DEBUG            "-O2 -g -Ob0 -noalign")
   set(CMAKE_C_FLAGS_DEBUGFULL        "-g -Ob0 -noalign")

#   set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -ansi -Wpointer-arith -fno-common")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -Wpointer-arith -fno-exceptions -fno-common")

   endif (CMAKE_C_COMPILER MATCHES "icc")
endif (CMAKE_SYSTEM_NAME MATCHES Linux)


#STRING(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)

SET (BLA CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER})

MESSAGE("--- Build Type: ${CMAKE_BUILD_TYPE}")
MESSAGE("--- Compiler Flags: ${CMAKE_CXX_FLAGS} ${${BLA}}")

ENDMACRO ( Check_Compiler )