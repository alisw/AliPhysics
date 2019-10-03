# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************

if(CAN_GENERATE_DOC)

  # Prepare Doxygen config from template
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    @ONLY
  )

  # Generate doc in build dir by using "make doxygen"
  add_custom_target(doxygen
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating AliPhysics reference manual using Doxygen"
  )

  # Install documentation in the installation directory if requested
  add_custom_target(install-doxygen
    rsync -a --delete "${CMAKE_CURRENT_BINARY_DIR}/html/" "${CMAKE_INSTALL_PREFIX}/doxygen/"
    DEPENDS doxygen
    COMMENT "Installing AliPhysics reference manual"
  )

endif(CAN_GENERATE_DOC)
