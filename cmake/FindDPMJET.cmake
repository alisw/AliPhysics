message(STATUS "looking for DPMJET in \"${DPMJET}"\")

find_library(DPMJET_LIB NAMES DPMJET PATHS ${DPMJET}/lib NO_DEFAULT_PATH)

if (DPMJET_LIB)
   set(DPMJET_LIBPATH ${DPMJET}/lib)
   set(DPMJET_FOUND TRUE)
  message(STATUS "DPMJET found")
else()
  message(STATUS "DPMJET not found")
endif()
