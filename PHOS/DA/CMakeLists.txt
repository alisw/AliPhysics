# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PHOS/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSbase
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSshuttle
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

set(STATIC_DEPENDENCIES PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHS" "BCM" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHS" "GAIN" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHS" "LED" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHS" "PED" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("CPV" "PED" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("CPV" "BCM" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("CPV" "GAIN" "${STATIC_DEPENDENCIES}")
