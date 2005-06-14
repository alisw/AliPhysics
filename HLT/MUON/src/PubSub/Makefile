
#BASEDIR := /home/HLT/src/testing
# Before Configure.py:
# BASEDIR := /home/HLT/src/current
BASEDIR := /home/dHLT/HLT-Data-Challenge-200411/src
ALIHLTDIR := $(BASEDIR)/AliHLT
#KIPDIR := $(BASEDIR)/kip
#MLUCDIR := $(KIPDIR)/MLUC
#BCLDIR := $(KIPDIR)/BCL
#PSIDIR := $(BASEDIR)/psi


# Before Configure.py:
# INTEGRTEST_6_ANALYSIS_DIR := /home/HLT/src/Util/Analysis/IntegrTest-6
INTEGRTEST_6_ANALYSIS_DIR := /home/dHLT/HLT-Data-Challenge-200411/src/IntegrTest-6

ARCH := $(shell uname -m)
SYSTEM := $(shell uname -s)
TARGET := $(SYSTEM)-$(ARCH)

TOPDIRREL := $(ALIHLTDIR)
MLUCTOPDIR := $(MLUCDIR)
BCLTOPDIR := $(BCLDIR)

include $(ALIHLTDIR)/conf/Makefile.conf

include $(MLUCTOPDIR)/conf/Makefile.extern.conf
include $(BCLTOPDIR)/conf/Makefile.extern.conf

include $(ALIHLTDIR)/Base/conf/Makefile.extern.conf
include $(ALIHLTDIR)/SC/conf/Makefile.extern.conf
include $(ALIHLTDIR)/Trigger/conf/Makefile.extern.conf
include $(ALIHLTDIR)/PubSub/conf/Makefile.extern.conf
include $(ALIHLTDIR)/DDL/conf/Makefile.extern.conf
include $(ALIHLTDIR)/WorkerComp/conf/Makefile.extern.conf


INTEGRTEST_6_ANALYSIS_INCLUDE := -I$(INTEGRTEST_6_ANALYSIS_DIR)
INTEGRTEST_6_ANALYSIS_LIB := $(INTEGRTEST_6_ANALYSIS_DIR)/lib/$(TARGET)/libIntegrTest-6-Analysis.a

INCLUDE := $(INCLUDE) $(INTEGRTEST_6_ANALYSIS_INCLUDE) -I../../src
LIBS := $(LIBS) $(INTEGRTEST_6_ANALYSIS_LIB)
DEFINES := -Dno_root


vpath %.cpp ../ ../Tracking/

all: TrackerComponent

TrackerComponent : TrackerComponent.o TrackerCallback.o MansoTracker.o Calculations.o $(INTEGRTEST_6_ANALYSIS_LIB)
	$(LD) $+ $(LDOPT) $(LIBS) -o $@

TrackerComponent.o : TrackerComponent.cpp TrackerComponent.hpp ../Cluster.hpp
	$(CC) -c $(INCLUDE) $(DEFINES) $(CCOPT) -o $@ $<

TrackerCallback.o : TrackerCallback.cpp TrackerCallback.hpp ../Cluster.hpp
	$(CC) -c $(INCLUDE) $(DEFINES) $(CCOPT) -o $@ $<

%.o : %.cpp
	$(CC) -c $(INCLUDE) $(DEFINES) $(CCOPT) -o $@ $<

clean:
	rm *~*
	rm -f *.o
	rm TrackerComponent
