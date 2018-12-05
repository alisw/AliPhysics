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
set(MODULE HIJING)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/${MODULE}
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(SYSTEM ${ROOT_INCLUDE_DIR})
include_directories(${AliRoot_SOURCE_DIR}/HIJING/hijing
                   )

# Sources
set(SRCS
  AliHijingRndm.cxx
)

set(FSRCS
hijing1_36/ar3jet.F hijing1_36/aran9.F hijing1_36/arorie.F hijing1_36/atrobo.F hijing1_36/attflv.F hijing1_36/attrad.F hijing1_36/bk.F hijing1_36/crsjet.F hijing1_36/fgp1.F hijing1_36/fgp2.F hijing1_36/fgp3.F hijing1_36/fhin.F hijing1_36/fjet.F hijing1_36/fjetrig.F hijing1_36/flap.F hijing1_36/flap1.F hijing1_36/flap2.F hijing1_36/fnjet.F hijing1_36/fnkick.F hijing1_36/fnkick2.F hijing1_36/fnstru.F hijing1_36/fnstrum.F hijing1_36/fnstrus.F hijing1_36/ftot.F hijing1_36/ftotjet.F hijing1_36/ftotrig.F hijing1_36/g.F hijing1_36/gauss1.F hijing1_36/gauss2.F hijing1_36/gauss3.F hijing1_36/gauss4.F hijing1_36/ghvq.F hijing1_36/gmin.F hijing1_36/gmre.F hijing1_36/gphoton.F hijing1_36/hiboost.F hijing1_36/hidata.F hijing1_36/hifun.F hijing1_36/hijcrs.F hijing1_36/hijcsc.F hijing1_36/hijels.F hijing1_36/hijflv.F hijing1_36/hijfrg.F hijing1_36/hijhrd.F hijing1_36/hijing.F hijing1_36/hijini.F hijing1_36/hijset.F hijing1_36/hijsft.F hijing1_36/hijsrt.F hijing1_36/hijwds.F hijing1_36/hiptdi.F hijing1_36/hirnd.F hijing1_36/hirnd2.F hijing1_36/jetini.F hijing1_36/omg0.F hijing1_36/parton.F hijing1_36/profile.F hijing1_36/quench.F hijing1_36/romg.F hijing1_36/rwdsax.F hijing1_36/sgmin.F hijing1_36/subcrs1.F hijing1_36/subcrs2.F hijing1_36/subcrs3.F hijing1_36/subcrs4.F hijing1_36/subcrs5.F hijing1_36/subcrs6.F hijing1_36/subcrs7.F hijing1_36/title.F hijing1_36/vegas.F hijing1_36/wdsax.F hijing1_36/wdsax1.F hijing1_36/wdsax2.F hipyset1_35/hirobo.F hipyset1_35/klu_hijing.F hipyset1_35/lu1ent_hijing.F hipyset1_35/lu2ent_hijing.F hipyset1_35/lu3ent_hijing.F hipyset1_35/lu4ent_hijing.F hipyset1_35/luboei_hijing.F hipyset1_35/lucell_hijing.F hipyset1_35/luchge_hijing.F hipyset1_35/luclus_hijing.F hipyset1_35/lucomp_hijing.F hipyset1_35/ludata_hijing.F hipyset1_35/ludecy_hijing.F hipyset1_35/luedit_hijing.F hipyset1_35/lueevt_hijing.F hipyset1_35/luerrm_hijing.F hipyset1_35/luexec_hijing.F hipyset1_35/lufowo_hijing.F hipyset1_35/lugive_hijing.F hipyset1_35/luhepc_hijing.F hipyset1_35/luindf_hijing.F hipyset1_35/lujmas_hijing.F hipyset1_35/lujoin_hijing.F hipyset1_35/lukfdi_hijing.F hipyset1_35/lulist_hijing.F hipyset1_35/luname_hijing.F hipyset1_35/luonia_hijing.F hipyset1_35/luprep_hijing.F hipyset1_35/luptdi_hijing.F hipyset1_35/luradk_hijing.F hipyset1_35/lurobo_hijing.F hipyset1_35/lushow_hijing.F hipyset1_35/lusphe_hijing.F hipyset1_35/lustrf_hijing.F hipyset1_35/lutabu_hijing.F hipyset1_35/lutest_hijing.F hipyset1_35/luthru_hijing.F hipyset1_35/luupda_hijing.F hipyset1_35/lux3jt_hijing.F hipyset1_35/lux4jt_hijing.F hipyset1_35/luxdif_hijing.F hipyset1_35/luxjet_hijing.F hipyset1_35/luxkfl_hijing.F hipyset1_35/luxtot_hijing.F hipyset1_35/luzdis_hijing.F hipyset1_35/plu_hijing.F hipyset1_35/pydata_hijing.F hipyset1_35/pydiff_hijing.F hipyset1_35/pyfram_hijing.F hipyset1_35/pygamm_hijing.F hipyset1_35/pyi3au_hijing.F hipyset1_35/pyinit_hijing.F hipyset1_35/pyinki_hijing.F hipyset1_35/pyinre_hijing.F hipyset1_35/pykcut_hijing.F hipyset1_35/pyklim_hijing.F hipyset1_35/pykmap_hijing.F hipyset1_35/pymaxi_hijing.F hipyset1_35/pymult_hijing.F hipyset1_35/pyovly_hijing.F hipyset1_35/pyrand_hijing.F hipyset1_35/pyremn_hijing.F hipyset1_35/pyresd_hijing.F hipyset1_35/pyscat_hijing.F hipyset1_35/pysigh_hijing.F hipyset1_35/pyspen_hijing.F hipyset1_35/pyspli_hijing.F hipyset1_35/pysspa_hijing.F hipyset1_35/pystat_hijing.F hipyset1_35/pystfe_hijing.F hipyset1_35/pystfu_hijing.F hipyset1_35/pytest_hijing.F hipyset1_35/pythia_hijing.F hipyset1_35/pyw1au_hijing.F hipyset1_35/pyw2au_hijing.F hipyset1_35/pywidt_hijing.F hipyset1_35/pyxtot_hijing.F hipyset1_35/ulalps_hijing.F hipyset1_35/ulangl_hijing.F hipyset1_35/ulmass_hijing.F
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
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx ${FSRCS})

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
