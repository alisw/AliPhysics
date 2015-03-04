# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PMD/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/PMD/PMDrec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )


setDAflags()

set(STATIC_DEPENDENCIES CDB-static PMDrec-static RAWDatarec-static STEERBase-static)
generateDA("PMD" "GAIN" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES CDB-static PMDrec-static RAWDatarec-static STEERBase-static)
generateDA("PMD" "PEDESTAL" "${STATIC_DEPENDENCIES}")
