############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

BINARY  = tracker

SOURCES = tracker_main.cxx \
	System/Mutex.cxx \
	System/SystemError.cxx \
	Version.cxx \
	Error.cxx \
	Utils.cxx \
	RegionOfInterest.cxx \
	Debug/ClusterSource.cxx \
	Debug/print.cxx \
	Tracking/IOHandler.cxx \
	Tracking/EventHandler.cxx

MACROS = 

LIBRARY_PATHS = 
LIBRARIES = dHLT

############################################################################

include build/config.mk
include $(BUILD_DIR)/rules.mk
