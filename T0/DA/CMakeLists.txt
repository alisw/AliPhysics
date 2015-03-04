# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/T0/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/T0/T0base
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

set(STATIC_DEPENDENCIES RAWDatarec-static T0base-static)
generateDA("T00" "Cosmic" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static T0base-static)
generateDA("T00" "Laser" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES CDB-static RAWDatarec-static T0base-static)
generateDA("T00" "Phys" "${STATIC_DEPENDENCIES}")
