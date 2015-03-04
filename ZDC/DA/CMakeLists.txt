# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/ZDC/DA)

# Additional include folders in alphabetical order
include_directories(
                   )

set(STATIC_DEPENDENCIES RAWDatabase-static RAWDatarec-static ZDCbase-static)
generateDA("ZDC" "CALIBMB" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatabase-static RAWDatarec-static ZDCbase-static)
generateDA("ZDC" "EMD" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES  RAWDatabase-static RAWDatarec-static ZDCbase-static)
generateDA("ZDC" "LASER" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatabase-static RAWDatarec-static ZDCbase-static)
generateDA("ZDC" "MAPPING" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatabase-static RAWDatarec-static ZDCbase-static)
generateDA("ZDC" "PEDESTAL" "${STATIC_DEPENDENCIES}")

