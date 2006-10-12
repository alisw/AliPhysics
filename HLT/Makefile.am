# $Id$
# Makefile template for the Alice HLT framework and components


if EN_HLT_SAMPLE
SAMPLE_DIR=SampleLib	
endif

if EN_HLT_TPC
TPC_DIR=TPCLib	
endif

if EN_HLT_PHOS
PHOS_DIR=PHOS
endif

SUBDIRS 		= BASE \
			  $(SAMPLE_DIR) \
			  $(TPC_DIR) \
			  $(PHOS_DIR) \
			  doc

DIST_SUBDIRS		= $(SUBDIRS)

#
# EOF
#
