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

#pragma link C++ class  AliITS+;
#pragma link C++ class  AliITSvPPRcoarseasymm+;
#pragma link C++ class  AliITSvPPRasymmFMD+;
#pragma link C++ class  AliITSvBeamTestITS04+;
#pragma link C++ class  AliITSvSPD02+;
#pragma link C++ class  AliITSvSDD03+;
#pragma link C++ class  AliITSvtest+;
#pragma link C++ class  AliITSvSSD03+;
#pragma link C++ class  AliITSv11+;
#pragma link C++ class  AliITSv11GeomCable+;
#pragma link C++ class  AliITSv11GeomCableFlat+;
#pragma link C++ class  AliITSv11GeomCableRound+;
#pragma link C++ class  AliITSv11Geometry+;
//#pragma link C++ class  AliITSv11GeometrySupport+;
//#pragma link C++ class  AliITSv11GeometrySPD+;
#pragma link C++ class  AliITSv11GeometrySDD+;
//#pragma link C++ class  AliITSv11GeometrySSD+;
#pragma link C++ class  AliITShit+;
// Standard ITS detector class initilizers
#pragma link C++ class  AliITSmodule+;
#pragma link C++ class  AliITSsimulation+;
#pragma link C++ class  AliITSsimulationSPD+;
#pragma link C++ class  AliITSsimulationSDD+;
#pragma link C++ class  AliITSsimulationSSD+;
#pragma link C++ class  AliITSTableSSD+;
#pragma link C++ class  AliITSsimulationFastPoints+;
#pragma link C++ class  AliITSsimulationFastPointsV0+;
#pragma link C++ class  AliITSDetTypeSim+;
#pragma link C++ class  AliITSstatistics+;
#pragma link C++ class  AliITSstatistics2+;
// These streamers must be formatted according to the raw data fromat

//        #pragma link C++ class  AliITSHNode+;
#pragma link C++ class  AliITSHuffman+;
#pragma link C++ class  AliITSetfSDD+;
// SSD simulation and reconstruction
#pragma link C++ class  AliITSdcsSSD+;
#pragma link C++ class AliITSsDigitize+;
#pragma link C++ class AliITSDigitizer+;
#pragma link C++ class AliITSFDigitizer+;
//#pragma link C++ class DisplayITSv11+;
// Raw data

#pragma link C++ class AliITSTrigger+;



#endif
