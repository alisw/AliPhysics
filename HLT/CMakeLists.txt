# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/EMCAL
${CMAKE_SOURCE_DIR}/HLT
${CMAKE_SOURCE_DIR}/HLT/BASE
${CMAKE_SOURCE_DIR}/HLT/BASE/HOMER
${CMAKE_SOURCE_DIR}/HLT/BASE/interface
${CMAKE_SOURCE_DIR}/HLT/BASE/util
${CMAKE_SOURCE_DIR}/HLT/EMCAL
${CMAKE_SOURCE_DIR}/HLT/ITS
${CMAKE_SOURCE_DIR}/HLT/MUON
${CMAKE_SOURCE_DIR}/HLT/MUON/OfflineInterface
${CMAKE_SOURCE_DIR}/HLT/MUON/OnlineAnalysis
${CMAKE_SOURCE_DIR}/HLT/MUON/utils
${CMAKE_SOURCE_DIR}/HLT/PHOS
${CMAKE_SOURCE_DIR}/HLT/RCU
${CMAKE_SOURCE_DIR}/HLT/SampleLib
${CMAKE_SOURCE_DIR}/HLT/TPCLib
${CMAKE_SOURCE_DIR}/HLT/TPCLib/comp
${CMAKE_SOURCE_DIR}/HLT/TPCLib/offline
${CMAKE_SOURCE_DIR}/HLT/TPCLib/tracking-ca
${CMAKE_SOURCE_DIR}/HLT/TRD
${CMAKE_SOURCE_DIR}/HLT/comp
${CMAKE_SOURCE_DIR}/HLT/offline
${CMAKE_SOURCE_DIR}/HLT/rec
${CMAKE_SOURCE_DIR}/HLT/shuttle
${CMAKE_SOURCE_DIR}/HLT/sim
${CMAKE_SOURCE_DIR}/HLT/tracking-ca
${CMAKE_SOURCE_DIR}/HLT/trigger
${CMAKE_SOURCE_DIR}/ITS
${CMAKE_SOURCE_DIR}/MUON
${CMAKE_SOURCE_DIR}/MUON/mapping
${CMAKE_SOURCE_DIR}/PHOS
${CMAKE_SOURCE_DIR}/RAW
${CMAKE_SOURCE_DIR}/STEER
${CMAKE_SOURCE_DIR}/TPC
${CMAKE_SOURCE_DIR}/TRD
${CMAKE_SOURCE_DIR}/global
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

ADD_DEFINITIONS(-Duse_aliroot -Duse_root -DUSE_ROOT -DROWHOUGHPARAMS -Duse_reconstruction -Duse_newio -DROOTVERSION=\"unchecked\" -DALIROOTVERSION=\"unchecked\" -D__ROOT__ -DUSE_ALILOG -DLINUX -DLIBHLTSIM_VERSION=0)

SetModule()

Include(CMake_libAliHLTHOMER.txt)

Include(CMake_libHLTinterface.txt)

Include(CMake_libAliHLTUtil.txt)

Include(CMake_libAliHLTSample.txt)

Include(CMake_libAliHLTRCU.txt)

Include(CMake_libHLTshuttle.txt)

Include(CMake_libAliHLTTrigger.txt)

Include(CMake_libAliHLTMUON.txt)

Include(CMake_libHLTrec.txt)

Include(CMake_libHLTbase.txt)

Include(CMake_libAliHLTGlobal.txt)

Include(CMake_libAliHLTTPC.txt)

Include(CMake_libAliHLTTRD.txt)

Include(CMake_libAliHLTComp.txt)

Include(CMake_libAliHLTITS.txt)

Include(CMake_libHLTsim.txt)

Include(CMake_libAliHLTPHOS.txt)

Include(CMake_libAliHLTEMCAL.txt)

Include(CMake_bindHLTdumpraw.txt)

