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

# Module
set(MODULE HERWIG)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/${MODULE})

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR})

# Sources
set(SRCS
  AliHerwigRndm.cxx
)

set (CSRCS  main.c herwig6_address.c)

set(FSRCS
  herwig6_called_from_cc.f herwig6_common_block_address.f herwig6510.f hwlhin.f str.f jimmy/cernmisc/jmabend.F jimmy/cernmisc/jmdbs4.F jimmy/cernmisc/jmdbsir3.F jimmy/cernmisc/jmdbska.F jimmy/cernmisc/jmlnoc.F jimmy/cernmisc/jmmlpt.F jimmy/cernmisc/jmmlst.F jimmy/cernmisc/jmucpy.F jimmy/divon4/addbnd.F jimmy/divon4/bfgs.F jimmy/divon4/bndopt.F jimmy/divon4/bndtst.F jimmy/divon4/bounds.F jimmy/divon4/bucmve.F jimmy/divon4/bufopt.F jimmy/divon4/delbnd.F jimmy/divon4/delete.F jimmy/divon4/delslv.F jimmy/divon4/divon.F jimmy/divon4/dvcopy.F jimmy/divon4/dvdot.F jimmy/divon4/dvnbkd.F jimmy/divon4/exmbuc.F jimmy/divon4/feasmv.F jimmy/divon4/feqn.F jimmy/divon4/fun.F jimmy/divon4/genpnt.F jimmy/divon4/grdcmp.F jimmy/divon4/intgrl.F jimmy/divon4/ldlsol.F jimmy/divon4/locsch.F jimmy/divon4/modchl.F jimmy/divon4/mulchk.F jimmy/divon4/mxstep.F jimmy/divon4/newptq.F jimmy/divon4/nmdchl.F jimmy/divon4/nocut.F jimmy/divon4/nodaud.F jimmy/divon4/orthvc.F jimmy/divon4/partn.F jimmy/divon4/quad.F jimmy/divon4/quasi.F jimmy/divon4/ranums.F jimmy/divon4/recparz.F jimmy/divon4/rlen.F jimmy/divon4/settol.F jimmy/divon4/shrnk.F jimmy/divon4/split.F jimmy/divon4/treaud.F jimmy/divon4/tredmp.F jimmy/divon4/tstext.F jimmy/divon4/usrint.F jimmy/divon4/usrtrm.F jimmy/jimmy/dfun.F jimmy/jimmy/dvnopt.F jimmy/jimmy/eikfun.F jimmy/jimmy/eikon.F jimmy/jimmy/epfun.F jimmy/jimmy/hwmsct.F jimmy/jimmy/hwnsct.F jimmy/jimmy/jimmin.F jimmy/jimmy/jmarea.F jimmy/jimmy/jmefin.F jimmy/jimmy/jmfact.F jimmy/jimmy/jmgamz.F jimmy/jimmy/jmgrid.F jimmy/jimmy/jminit.F jimmy/jimmy/jmkern.F jimmy/jimmy/jmknif.F jimmy/jimmy/jmpn.F jimmy/jimmy/jmrndm.F jimmy/jimmy/jmsn.F jimmy/jimmy/jmsnfn.F jimmy/jimmy/jmxs1.F jimmy/jimmy/qcdres.F
)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS MathCore)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx ${FSRCS} ${CSRCS})

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

message(STATUS "${MODULE} enabled")
