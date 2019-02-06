# Include the libraries
 
# Fortran modules
if(CMAKE_Fortran_COMPILER)
  add_subdirectory(AliFemto)
endif(CMAKE_Fortran_COMPILER)

add_subdirectory(AliFemtoUser)
add_subdirectory(Chaoticity)
add_subdirectory(ESE)
add_subdirectory(HBTwrtEP)
add_subdirectory(HCascadeAnalysis)
add_subdirectory(K0Analysis)
add_subdirectory(UNICOR)
add_subdirectory(V0LamAnalysis)
add_subdirectory(PLamAnalysisPP)
add_subdirectory(FemtoDream)
add_subdirectory(KPAnalysis)
add_subdirectory(EfficiencyTask)

install(DIRECTORY macros DESTINATION PWGCF/FEMTOSCOPY)
