# For more instructions please read INSTALLDA.txt

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/AD/DA)

# Additional include folders in alphabetical order
include_directories(
		)
		   
set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static ADbase-static)
# Generate the DA using DetectorName, Algortihm
# The DA file has to be named DetectorNameAlgortihmda.cxx
# generateDA("DetectorName" "Algorithm" "${STATIC_DEPENDENCIES}")
generateDA("AD0" "PEDESTAL" "${STATIC_DEPENDENCIES}")
generateDA("AD0" "TIMESLEWING" "${STATIC_DEPENDENCIES}")
