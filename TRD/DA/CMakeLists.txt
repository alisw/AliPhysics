# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TRD/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCbase
                   )

set(STATIC_DEPENDENCIES CDB-static RAWDatarec-static STEERBase-static TRDbase-static)
generateDA("TRD" "HALFCHAMBERSTATUS" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES CDB-static RAWDatarec-static STEERBase-static TRDbase-static)
generateDA("TRD" "PEDESTAL" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TRDbase-static)
generateDA("TRD" "VDRIFT" "${STATIC_DEPENDENCIES}")
