
Flow ReadMe

The Flow package in AliRoot consistes of two parts.

1) Macros for fast flow analysis used for PPR studies
2) Interface to STAR flow analysis package


>> Alice PPR fast flow analysis package

This package is composed out of macros:

AliFlowPicoCreate.C - create pico event structures in "flowPicoEvent.root"
AliFlowReconstruction.C - fill pico event structure from Kinematics

All analysis macros are using "flowPicoEvent.root" file with data

AliFlowDist.C  - draw different distributions 
		 for a given multiplicity and V2

AliFlowDrawV2.C - Draw event plane resolution and V2 resolution
		  for a given multiplicity in function of V2

AliFlowDrawSummary.C - Draw resolution in function of multiplicity 
		       and V2.

>> Alice Interface to STAR Flow analysis package


