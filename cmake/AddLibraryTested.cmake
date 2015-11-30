function(add_library_tested NAME ...)
  # add library
  add_library(${ARGV})

  # and remember it for the list of libraries
  list(APPEND ALILIBSTESTED ${NAME})
  set(ALILIBSTESTED ${ALILIBSTESTED} CACHE INTERNAL "ALILIBSTESTED")
endfunction(add_library_tested)
