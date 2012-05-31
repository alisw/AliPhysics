#/**************************************************************************
# * This file is property of and copyright by the ALICE HLT Project        * 
# * All rights reserved.                                                   *
# *                                                                        *
# * Primary Authors:                                                       *
# *   Artur Szostak <artursz@iafrica.com>                                  *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          * 
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************/

# This makefile allows the dHLTdumpraw utility to be built standalone without
# requiring external libraries. All system libraries are built with -static.
# It may be useful in certain cases where no AliRoot libraries are available or
# it is not certain which library versions are in use.
# In general and usual cases one should rather use the utilities as built by
# the normal AliRoot build system since this makefile is less portable between
# platforms.

SRCS =	AliHLTMUONUtils.cxx \
	AliHLTMUONConstants.cxx \
	AliHLTMUONTriggerRecordsBlockStruct.cxx \
	AliHLTMUONTrigRecsDebugBlockStruct.cxx \
	AliHLTMUONTriggerChannelsBlockStruct.cxx \
	AliHLTMUONRecHitsBlockStruct.cxx \
	AliHLTMUONClustersBlockStruct.cxx \
	AliHLTMUONChannelsBlockStruct.cxx \
	AliHLTMUONMansoTracksBlockStruct.cxx \
	AliHLTMUONMansoCandidatesBlockStruct.cxx \
	AliHLTMUONSinglesDecisionBlockStruct.cxx \
	AliHLTMUONPairsDecisionBlockStruct.cxx

vpath %.h ../
vpath %.cxx ../

HDRS := $(SRCS:.cxx=.h)
DEPS := $(SRCS:.cxx=.d)
OBJS := $(SRCS:.cxx=.o)

.PHONY : all clean

all : dHLTdumpraw

clean :
	rm -f *.o *.d

dHLTdumpraw : $(OBJS) dHLTdumpraw.o
	g++ -O -static -static-libgcc $^ -o $@

%.o : %.cxx
	g++ -O -I../../BASE -I../ -Wall -Weffc++ -c $< -o $@

SEDCMD := "s/\(.*\)\.o/$(subst /,\/,$(@D)/)\1\.o $(subst /,\/,$(@D)/)\1\.d/g"

%.d : %.cxx
	g++ -I../../BASE -I../ -MM $< | sed $(SEDCMD) > $@

-include $(DEPS)
