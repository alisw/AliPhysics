# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MUON/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcalib
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcore
                    ${AliRoot_SOURCE_DIR}/MUON/MUONmapping
                    ${AliRoot_SOURCE_DIR}/MUON/MUONraw
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )


set(STATIC_DEPENDENCIES MUONcalib-static MUONcore-static MUONmapping-static MUONraw-static RAWDatarec-static)
generateDA("MTR" "" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES MUONcalib-static MUONcore-static MUONmapping-static MUONraw-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("MCH" "OCC" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES MUONcalib-static MUONcore-static MUONmapping-static MUONraw-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("MCH" "BPEVO" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ESD-static MUONcalib-static MUONcore-static MUONmapping-static MUONraw-static RAWDatarec-static STEERBase-static)
generateDA("MCH" "PED" "${STATIC_DEPENDENCIES}")
