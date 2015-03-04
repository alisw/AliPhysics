# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/ACORDE/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/ACORDE/ACORDEbase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/dateStream
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                   )

set(STATIC_DEPENDENCIES  ESD-static RAWDatabase-static RAWDatarec-static ACORDEbase-static)
generateDA("ACO" "" "${STATIC_DEPENDENCIES}")