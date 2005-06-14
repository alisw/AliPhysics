############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

BINARY  = signalL2

SOURCES = signalL2.cxx \
	System/SystemError.cxx \
	System/Socket.cxx \
	System/Mutex.cxx \
	System/MutexCondition.cxx \
	System/Thread.cxx \
	System/SignalHandler.cxx \
	BCMP/Packets.cxx \
	BCMP/Sender.cxx \
	BCMP/EventQueue.cxx \
	DDL/L2SignalSender.cxx \
	Debug/print.cxx \
	Version.cxx \
	Utils.cxx \
	new.cxx \
	Error.cxx

MACROS =

LIBRARY_PATHS = 
LIBRARIES = pthread

############################################################################

include build/config.mk
include $(BUILD_DIR)/rules.mk
