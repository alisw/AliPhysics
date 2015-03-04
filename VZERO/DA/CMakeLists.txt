# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/VZERO/DA)

# Additional include folders in alphabetical order
include_directories(
                   )


set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static VZEROrec-static)
generateDA("V00" "CHANNELDELAYTUNING" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static VZEROrec-static)
generateDA("V00" "" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static VZEROrec-static)
generateDA("V00" "INTEGRATIONGATETUNING" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static VZEROrec-static)
generateDA("V00" "PbPb" "${STATIC_DEPENDENCIES}")
