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
#pragma link C++ class  AliITSv11+;
#pragma link C++ class  AliITSv11GeomCable+;
#pragma link C++ class  AliITSv11GeomCableFlat+;
#pragma link C++ class  AliITSv11GeomCableRound+;
#pragma link C++ class  AliITSv11Geometry+;
#pragma link C++ class  AliITSv11GeometrySupport+;
#pragma link C++ class  AliITSv11GeometrySPD+;
#pragma link C++ class  AliITSv11GeometrySDD+;
#pragma link C++ class  AliITSv11GeometrySSD+;
#pragma link C++ class  AliITSv11Hybrid+;
#pragma link C++ class  AliITShit+;
// Standard ITS detector class initilizers
#pragma link C++ class  AliITSmodule+;
#pragma link C++ class  AliITSsimulation+;
#pragma link C++ class  AliITSsimulationSPD+;
#pragma link C++ class  AliITSsimulationSDD+;
#pragma link C++ class  AliITSsimulationSSD+;
#pragma link C++ class  AliITSTableSSD+;
#pragma link C++ class  AliITSsimulationFastPoints+;
#pragma link C++ class  AliITSSimuParam+;
#pragma link C++ class  AliITSDetTypeSim+;

#pragma link C++ class  AliITSetfSDD+;
// SSD simulation and reconstruction
#pragma link C++ class AliITSsDigitize+;
#pragma link C++ class AliITSDigitizer+;
//#pragma link C++ class DisplayITSv11+;
// Raw data


#pragma link C++ class AliITSFOEfficiencySPD+;
#pragma link C++ class AliITSFOEfficiencySPDColumn+;
#pragma link C++ class AliITSFONoiseSPD+;
#pragma link C++ class AliITSFOGeneratorSPD+;
#pragma link C++ class AliITSFOSignalsSPD+;

#pragma link C++ class AliITSTrigger+;
#pragma link C++ class AliITSTriggerFOProcessor+;
#pragma link C++ class AliITSQADataMakerSim+;
#pragma link C++ class AliITSQASPDDataMakerSim+;
#pragma link C++ class AliITSQASDDDataMakerSim+;
#pragma link C++ class AliITSQASSDDataMakerSim+;




#endif
