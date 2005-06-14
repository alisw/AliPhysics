############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

LIBRARY = dHLT

SOURCES = Version.cxx

INCLUDE_HEADERS = Version.hpp BasicTypes.hpp

MACROS =

#########################################################################

include build/config.mk
include $(BUILD_DIR)/rules.mk

