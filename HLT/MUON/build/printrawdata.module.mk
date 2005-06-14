############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

BINARY  = printrawdata

SOURCES = Debug/printrawdata.cxx \
	System/File.cxx \
	System/SystemError.cxx \
	RegionOfInterest.cxx \
	Version.cxx \
	Utils.cxx \
	new.cxx \
	Error.cxx

MACROS =

############################################################################

include build/config.mk
include $(BUILD_DIR)/rules.mk
