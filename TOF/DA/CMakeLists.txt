# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TOF/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/TOF/TOFbase
                   )

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TOFbase-static)
generateDA("TOF" "calib" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TOFbase-static)
generateDA("TOF" "" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TOFbase-static)
generateDA("TOF" "noise" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TOFbase-static)
generateDA("TOF" "physics" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TOFbase-static)
generateDA("TOF" "pulser" "${STATIC_DEPENDENCIES}")

