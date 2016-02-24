/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file MUONsimLinkDef.h
/// \brief The CINT link definitions for \ref sim 

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ class AliMUON+; 
#pragma link C++ class AliMUONv1+; 
#pragma link C++ class AliMUONHit+; 

// geometry 
#pragma link C++ class AliMUONCommonGeometryBuilder+;
#pragma link C++ class AliMUONSt1SpecialMotif+; 
#pragma link C++ class AliMUONSt1GeometryBuilder+; 
#pragma link C++ class AliMUONSt1GeometryBuilderV2+; 
#pragma link C++ class AliMUONSt2GeometryBuilder+; 
#pragma link C++ class AliMUONSt2GeometryBuilderV2+; 
#pragma link C++ class AliMUONSlatGeometryBuilder+; 
#pragma link C++ class AliMUONTriggerGeometryBuilder+; 

// digitizer
#pragma link C++ class AliMUONDigitizerV3+; 
#pragma link C++ class AliMUONSDigitizerV2+;  
#pragma link C++ class AliMUONVHitStore+;
#pragma link C++ class AliMUONHitStoreV1+;

// response
#pragma link C++ class AliMUONChamber+; 
#pragma link C++ class AliMUONResponse+; 
#pragma link C++ class AliMUONResponseV0+;
#pragma link C++ class AliMUONResponseTrigger+; 
#pragma link C++ class AliMUONResponseTriggerV1+;
#pragma link C++ class AliMUONResponseFactory+;  

// trigger
#pragma link C++ class AliMUONTrigger+;
#pragma link C++ class AliMUONChamberTrigger+; 

// misc
#pragma link C++ class AliMUONMCDataInterface+;
#pragma link C++ class AliMUONPedestalEventGenerator+;

#pragma link C++ class AliMUONQADataMakerSim+;

#endif


