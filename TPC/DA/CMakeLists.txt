# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TPC/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCbase
                   )

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TPCbase-static)
generateDA("TPC" "CE" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TPCbase-static)
generateDA("TPC" "CEnew" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TPCbase-static)
generateDA("TPC" "PEDESTAL" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TPCbase-static)
generateDA("TPC" "PULSER" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TPCbase-static)
generateDA("TPC" "QA" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TPCbase-static)
generateDA("TPC" "RAW" "${STATIC_DEPENDENCIES}")
