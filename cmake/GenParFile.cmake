macro(add_parfile PARMODULE PARSOURCES PARHEADERS PARLINKDEF PARLIBDEPS)

  message(STATUS "PARfile generation: ${PARMODULE}")

  # Libraries
  foreach(THISLIB ${PARLIBDEPS})
    set(_PARLIBDEPS "${_PARLIBDEPS} lib${THISLIB}")
  endforeach()
  string(STRIP "${_PARLIBDEPS}" PARLIBDEPS)

  # Export variable
  set(PARMODULE "${PARMODULE}")

  message(STATUS "--> Module: ${PARMODULE}")
  message(STATUS "--> Sources: ${PARSOURCES}")
  message(STATUS "--> Deps: ${PARLIBDEPS}")

  foreach(LOOPVAR ${PARSOURCES})
    message(STATUS "----> ${LOOPVAR}")
  endforeach()

  # PARfile output directory (the one we will tar)
  set(PARDIR ${CMAKE_CURRENT_BINARY_DIR}/PARfiles/${PARMODULE})

  # Create base directory for this module's PARfile: this is the directory we will tar
  # This works as "mkdir -p" (i.e. it's recursive and creates parents)
  file(MAKE_DIRECTORY ${PARDIR}/PROOF-INF)

  # Create Makefile
  configure_file(
      ${PROJECT_SOURCE_DIR}/cmake/PARfiles/Makefile.in
      ${PARDIR}/Makefile
      @ONLY
  )

  # Create BUILD.sh
  configure_file(
      ${PROJECT_SOURCE_DIR}/cmake/PARfiles/BUILD.sh.in
      ${PARDIR}/PROOF-INF/BUILD.sh
      @ONLY
  )
  execute_process(COMMAND chmod a+x ${PARDIR}/PROOF-INF/BUILD.sh)

  # Create SETUP.C
  configure_file(
      ${PROJECT_SOURCE_DIR}/cmake/PARfiles/SETUP.C.in
      ${PARDIR}/PROOF-INF/SETUP.C
      @ONLY
  )

  # Target for creating PARfile (would stop after the first failed COMMAND)
  add_custom_target("${PARMODULE}.par"
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && cp ${PARSOURCES} ${PARHEADERS} ${PARLINKDEF} ${PARDIR}/
    COMMAND cmake -E copy ${ROOT_ETCDIR}/Makefile.arch ${PARDIR}/
    COMMAND tar -C ${PARDIR}/.. -czf ${PARDIR}/../${PARMODULE}.par ${PARMODULE}/
  )

endmacro()
