############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

LIBRARY = MUONHLT

DICTIONARY_HEADERS = AliRoot/Region.hpp \
	AliRoot/Point.hpp \
	AliRoot/TriggerRecord.hpp \
	AliRoot/ADCStream.hpp \
	AliRoot/Track.hpp \
	AliRoot/ADCStreamSource.hpp \
	AliRoot/TriggerSource.hpp \
	AliRoot/ClusterSource.hpp \
	AliRoot/TrackSink.hpp \
	AliRoot/TrackerCallback.hpp \
	AliRoot/ClusterFinderCallback.hpp \
	AliRoot/MicrodHLT.hpp

STUB_HEADERS = AliRoot/TrackerInterface.hpp \
	AliRoot/ClusterFinderInterface.hpp

SOURCES = $(DICTIONARY_HEADERS:.hpp=.cxx)

SOURCES += Utils.cxx \
	Error.cxx \
	Version.cxx \
	RegionOfInterest.cxx \
	Tracking/Calculations.cxx \
	Tracking/MansoTracker.cxx \
	Clustering/CenterOfGravityFinder.cxx \
	AliRoot/TrackerProxy.cxx \
	AliRoot/ClusterFinderProxy.cxx \
	AliRoot/convert.cxx

ifdef DEBUGGING
SOURCES += Debug/print.cxx
endif

INCLUDE_HEADERS = $(DICTIONARY_HEADERS) AliRoot/Tracker.hpp

DICTIONARY  = Dictionary.cxx
LINKDEF     = AliRoot/MUONHLTLinkDef.hpp

INCLUDES = $(ALICE_ROOT)/MUON $(ALICE_ROOT)/include/ $(ROOTSYS)/include/
MACROS = __ROOT__

# Test to see if the AliLog.h file exists in the AliRoot distribution.
# If it does then prompt the code to use it.
ifneq ($(findstring AliLog.h,$(shell find $$ALICE_ROOT/include -name AliLog.h)),)
MACROS += USE_ALILOG
endif

ifndef DEBUGGING
MACROS += LOG_NO_DEBUG
endif

#########################################################################

include build/config.mk
include $(BUILD_DIR)/rules.mk
