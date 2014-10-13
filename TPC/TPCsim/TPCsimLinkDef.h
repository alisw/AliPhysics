#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliTPChit+;    // defined in AliTPC.h
                                      // Container: Space-point coordinates + charge + track label + VolumeID
                                      // --- move to extra class 

#pragma link C++ class AliTrackHitsParamV2+; // defined in AliTPCTrackHitsV2.h 
                                             // Local parabolic parametrization of hits in global coordinates
                                             // --- move to extra class 
#pragma link C++ class AliTPCTrackHitsV2+;   // Containter: For Tracks Hits ("compressed" AliTPChit)

// --- MOVE AliTPC -> AliTPCMC
// AliTPCv0 -> AliTPCMCCoarse
// AliTPCv1 -> AliTPCMCParametrized
// AliTPCv2 -> AliTPCMCDefault
// AliTPCv4 -> AliTPCMCKrypton
// AliTPCLaser -> AliTPCMCLaser

#pragma link C++ class AliTPC+;       // Base class: for TPC simulation: eg defines Material, Geometry, etc ...   
                                      // --- Check to move relevant parameters to AliTPCParam
#pragma link C++ class AliTPCv0+;     // Coarse geometry (no sensitive volume)
                                      // --- Update Documentation
#pragma link C++ class AliTPCv2+;     // Default version - is used
                                      // --- Update Documentation
#pragma link C++ class AliTPCv4+;     // Krypton simulation - is used
                                      // --- Update Documentation
#pragma link C++ class AliTPCLaser+;  // Laser Simulation 
                                      // --- Update Documentation

#pragma link C++ class AliTPCDigitizer;   // Create Digits out of SDigits and partially SDigits from Hits
                                          // --- Update Documentation

#pragma link C++ class AliTPCBuffer+;     // Used to to write digit in raw format
                                          // --- Update Documentation
#pragma link C++ class AliTPCDDLRawData+; // Used to to write digit in raw format
                                          // --- Update Documentation

#pragma link C++ class AliTPCQADataMakerSim+; // Offline QA - to be dropped ...

#endif

