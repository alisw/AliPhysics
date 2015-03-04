# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/FMD/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/FMD/
                    ${AliRoot_SOURCE_DIR}/FMD/FMDbase
                    ${AliRoot_SOURCE_DIR}/FMD/FMDutil
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

set(STATIC_DEPENDENCIES FMDbase-static FMDutil-static RAWDatarec-static STEERBase-static)
generateDA("FMD" "Base" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES FMDbase-static FMDutil-static)
generateDA("FMD" "Gain" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES  FMDbase-static FMDutil-static)
generateDA("FMD" "Pedestal" "${STATIC_DEPENDENCIES}")
