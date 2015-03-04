# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/EMCAL/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALbase
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALUtils
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )


set(STATIC_DEPENDENCIES EMCALbase-static RAWDatarec-static RAWDatabase-static STEERBase-static)
generateDA("EMC" "PED" "${STATIC_DEPENDENCIES}")
generateDA("EMC" "LED" "${STATIC_DEPENDENCIES}")
