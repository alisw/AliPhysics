############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

BINARY  = triggerDDL

SOURCES = triggerDDL.cxx \
	System/SystemError.cxx \
	System/SignalHandler.cxx \
	System/Socket.cxx \
	System/Mutex.cxx \
	System/MutexCondition.cxx \
	System/Thread.cxx \
	System/File.cxx \
	System/Directory.cxx \
	System/Routines.cxx \
	BCMP/Packets.cxx \
	BCMP/Receiver.cxx \
	BCMP/EventQueue.cxx \
	DDL/FileList.cxx \
	DDL/L2SignalReceiver.cxx \
	Debug/print.cxx \
	Version.cxx \
	Utils.cxx \
	new.cxx \
	Error.cxx

MACROS = USE_GETDENTS_SYSCALL

LIBRARY_PATHS = 
LIBRARIES = pthread

############################################################################

include build/config.mk
include $(BUILD_DIR)/rules.mk
