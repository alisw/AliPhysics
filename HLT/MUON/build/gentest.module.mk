############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

BINARY  = gentest

SOURCES = gentest_main.cxx \
	Version.cxx \
	Utils.cxx \
	Error.cxx \
	RegionOfInterest.cxx \
	Tracking/IOHandler.cxx \
	Tracking/EventHandler.cxx \
	Debug/ClusterSource.cxx \
	Debug/TriggerSource.cxx \
	Debug/DataGenerator.cxx \
	Debug/print.cxx \
	Framework/Global.cxx \
	Buffers/GarbageCollector.cxx \
	Tracking/MansoTracker.cxx \
	Tracking/Calculations.cxx

MACROS = 

#########################################################################

include build/config.mk
include $(BUILD_DIR)/rules.mk
