# Top level Makefile.inc for ALICE's AMORE
# Author: Filimon Roukoutakis

include $(shell root-config --prefix)/test/Makefile.arch

ROOT_CFLAGS	:= $(filter-out -I% -D%, $(shell root-config --cflags))
ROOT_CPPFLAGS	:= $(filter-out $(ROOT_CFLAGS), $(shell root-config --cflags))
ROOT_LDFLAGS	:= $(filter-out -l%, $(shell root-config --libs))
ROOT_LIBS	:= $(filter -l%, $(shell root-config --libs))
DATE_CFLAGS	:= 
DATE_RCPROXYCFLAGS	:= $(shell date-config --rcproxycflags)
DATE_CPPFLAGS	:= $(patsubst -Dlong64="long long",,$(shell date-config --cflags))
DATE_LDFLAGS	:= $(filter-out -l%,$(shell date-config --monitorlibs))
DATE_LIBS	:= $(filter -l%, $(shell date-config --monitorlibs))
DATE_RCPROXYLIBS	:= $(shell date-config --rcproxylibs)
#ALIROOT_CPPFLAGS:= -I$(ALICE_ROOT)/include -I${ALICE_ROOT}/RAW
#ALIROOT_CFLAGS	:= 
#ALIROOT_LDFLAGS	:= -L$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)
#ALIROOT_LIBS	:= -lSTEER -lESD -lRAWDatabase
MYSQL_CFLAGS	:= $(shell mysql_config --cflags)
MYSQL_CPPFLAGS	:= $(shell mysql_config --include)
MYSQL_LIBS	:= $(shell mysql_config --libs)
AMORE_CPPFLAGS 	:= $(shell amore-config --includes)
AMORE_CFLAGS	:= $(shell amore-config --cflags) #-pg #-pedantic
AMORE_LDFLAGS	:= $(shell amore-config --ldflags) #-pg
AMORE_LIBS	:= $(shell amore-config --libs)

CPP		:= g++ -E 
CPPFLAGS	:= $(ROOT_CPPFLAGS) $(DATE_CPPFLAGS) $(DATE_RCPROXYCFLAGS) $(ALIROOT_CPPFLAGS) $(MYSQL_CPPFLAGS) $(AMORE_CPPFLAGS)
CC		:= gcc -c 
CFLAGS		:= $(ROOT_CFLAGS) $(DATE_CFLAGS) $(ALIROOT_CFLAGS) $(MYSQL_CFLAGS) $(AMORE_CFLAGS)
CXX		:= g++ -c 
CXXFLAGS	:= $(CFLAGS)
LD		:= g++ 
LDFLAGS		:= $(SOFLAGS) $(ROOT_LDFLAGS) $(DATE_LDFLAGS) $(ALIROOT_LDFLAGS) $(AMORE_LDFLAGS)
LIBS		:= $(ROOT_LIBS) $(DATE_LIBS) $(DATE_RCPROXYLIBS) $(ALIROOT_LIBS) $(MYSQL_LIBS) $(AMORE_LIBS)

HdrSuf = h
#SrcSuf = cpp

#%.$(ObjSuf):%.$(SrcSuf)
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< 

#%.$(ObjSuf):%.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) $<
#%:%.$(ObjSuf)
#	$(LD) $(LDFLAGS) $(LIBS) $^

TOP_DIR = ${AMORE_SITE}
SRC_DIR = $(TOP_DIR)/src
BIN_DIR = $(TOP_DIR)/bin
LIB_DIR = $(TOP_DIR)/lib

AllSources := $(wildcard *.$(SrcSuf))
HdrSuf      = h
AllHeaders := $(wildcard *.$(HdrSuf))
AllObjects := $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(wildcard *.$(SrcSuf)))
