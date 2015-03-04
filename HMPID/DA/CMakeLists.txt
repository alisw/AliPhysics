# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/ACORDE/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/HMPID/HMPIDbase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarecOnline
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )


set(STATIC_DEPENDENCIES HMPIDbase-static RAWDatarecOnline-static RAWDatarec-static RAWDatasim-static)
generateDA("HMP" "" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES HMPIDbase-static RAWDatarecOnline-static RAWDatarec-static)
generateDA("HMP" "phys" "${STATIC_DEPENDENCIES}")
