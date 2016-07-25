# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HMPID/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/ITS/ITSbase
                    ${AliRoot_SOURCE_DIR}/ITS/ITSrec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )



set(STATIC_DEPENDENCIES ITSbase-static ITSrec-static RAWDatarec-static)
generateDA("SDD" "BAS" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static)
generateDA("SDD" "GAIN" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static)
generateDA("SDD" "INJ" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static)
generateDA("SPD" "FO" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static)
generateDA("SPD" "PHYS" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static)
generateDA("SPD" "SCAN" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static RAWDatarec-static STEERBase-static CDB-static)
generateDA("SPD" "VertexDiamond" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSbase-static)
generateDA("SSD" "PED" "${STATIC_DEPENDENCIES}")

