#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
//#pragma link C++ enum   Cluster_t;

//#pragma link C++ global gITSdisplay;  // global used by AliITSdisplay

// Standard ITS classes 
 
#pragma link C++ class  AliITS-;
#pragma link C++ class  AliITSv1-;
#pragma link C++ class  AliITSv3-;
#pragma link C++ class  AliITSv5-;
#pragma link C++ class  AliITShit-;
#pragma link C++ class  AliITSdigit+;
#pragma link C++ class  AliITSdigitSPD+;
#pragma link C++ class  AliITSdigitSDD+;
#pragma link C++ class  AliITSdigitSSD+;
#pragma link C++ class  AliITSTransientDigit+;
#pragma link C++ class  AliITSgeom-;
#pragma link C++ class  AliITSgeomSPD+;
#pragma link C++ class  AliITSgeomSDD+;
#pragma link C++ class  AliITSgeomSSD+;

#pragma link C++ class AliITSgeomSPD300-;
#pragma link C++ class AliITSgeomSPD425-;

#pragma link C++ class  AliITSmodule+;
#pragma link C++ class  AliITSRecPoint+;
#pragma link C++ class  AliITSRawCluster+;
#pragma link C++ class  AliITSRawClusterSPD+;
#pragma link C++ class  AliITSRawClusterSDD+;
#pragma link C++ class  AliITSRawClusterSSD+;
#pragma link C++ class  AliITSMap+;
#pragma link C++ class  AliITSMapA1+;
#pragma link C++ class  AliITSMapA2+;
#pragma link C++ class  AliITSsegmentation+;
#pragma link C++ class  AliITSsegmentationSPD+;
#pragma link C++ class  AliITSsegmentationSDD+;
#pragma link C++ class  AliITSsegmentationSSD+;
#pragma link C++ class  AliITSresponse+;
#pragma link C++ class  AliITSresponseSPD+;
#pragma link C++ class  AliITSresponseSDD+;
#pragma link C++ class  AliITSresponseSSD+;
#pragma link C++ class  AliITSsimulation+;
#pragma link C++ class  AliITSsimulationSPD+;
#pragma link C++ class  AliITSsimulationSDD+;
#pragma link C++ class  AliITSsimulationSSD+;
#pragma link C++ class  AliITSsimulationFastPoints+;
#pragma link C++ class  AliITSsimulationFastPointsV0+;
#pragma link C++ class  AliITSClusterFinder+;
#pragma link C++ class  AliITSClusterFinderSPD+;
#pragma link C++ class  AliITSClusterFinderSDD+;
#pragma link C++ class  AliITSClusterFinderSSD+;
#pragma link C++ class  AliITSDetType+;
#pragma link C++ class  AliITStrack+;
// SDD simulation
#pragma link C++ class  AliITSRawData+;
#pragma link C++ class  AliITSInStream-;
#pragma link C++ class  AliITSOutStream-;
#pragma link C++ class  AliITSHNode+;
#pragma link C++ class  AliITSHTable+;
#pragma link C++ class  AliITSetfSDD+;
// SSD simulation and reconstruction
#pragma link C++ class  AliITSdictSSD+;
#pragma link C++ class  AliITSdcsSSD+;
#pragma link C++ class  AliITSclusterSSD+;
#pragma link C++ class  AliITSpackageSSD+;
// New used for Alignment studdies
#pragma link C++ class  AliITSAlignmentTrack-;
#pragma link C++ class  AliITSAlignmentModule-;
#pragma link C function HitsTo;
#pragma link C function HitsToClustAl;
#pragma link C function FillGlobalPositions;
#pragma link C function PlotGeomChanges;
#pragma link C function FitAllTracks;
#pragma link C function FitVertexAll;
#pragma link C function OnlyOneGeometry;
#pragma link C function deleteClustAl;
#pragma link C++ class  AliITSstatistics-;
#pragma link C++ class  AliITSstatistics2-;
// New used for AliITSdisplay
//#pragma link C++ class  AliITSdisplay;
//#pragma link C++ class  AliITSDisplay;
//#pragma link C++ class  TInputDialog;   // MUST BE RENAMED
//#pragma link C function OpenFileDialog;
//#pragma link C function GetStringDialog;
//#pragma link C function GetIntegerDialog;
//#pragma link C function GetFloatDialog;
// New classes used for Tracking

// This class will always be for ITS only
#pragma link C++ class  AliITSvtest-;

#endif
