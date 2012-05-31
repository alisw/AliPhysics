/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file MUONbaseLinkDef.h
/// \brief The CINT link definitions for \ref base

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// geometry 
#pragma link C++ class AliMUONMathieson+; 

// info classes 
#pragma link C++ class AliMUONConstants+; 

// containers
#pragma link C++ class AliMUONVCluster+;
#pragma link C++ class AliMUONRawCluster+;
#pragma link C++ class AliMUONRawClusterV2+;
#pragma link C++ class AliMUONDigit+; 
#pragma link C++ class AliMUONVDigit+; 
#pragma link C++ class AliMUONRealDigit+; 
#pragma link C++ class AliMUONVDigitStore+;
#pragma link C++ class AliMUONDigitStoreV1+;
#pragma link C++ class AliMUONDigitStoreV1Iterator+;
#pragma link C++ class AliMUONDigitStoreVImpl+;
#pragma link C++ class AliMUONDigitStoreVImplIterator+;
#pragma link C++ class AliMUONDigitStoreV2R+;
#pragma link C++ class AliMUONDigitStoreV2S+;
#pragma link C++ class AliMUONTOTCAStoreIterator+;
#pragma link C++ class AliMUONTriggerCircuit+;
#pragma link C++ class AliMUONVTriggerStore+;
#pragma link C++ class AliMUONTriggerStoreV1+;
#pragma link C++ class AliMUONTriggerElectronics+;

// raw data
#pragma link C++ class AliMUONDigitMaker+;
#pragma link C++ class AliMUONRawWriter+;

// calibration
#pragma link C++ class AliMUONDigitCalibrator+;
#pragma link C++ class AliMUONPadStatusMapMaker+;

// needed by qa (sim and rec)
#pragma link C++ class AliMUONRecoParam+;
#pragma link C++ class AliMUONPadStatusMaker+;
#pragma link C++ class AliMUONVQAChecker+;
#pragma link C++ class AliMUONQAChecker+;
#pragma link C++ class AliMUONTrackerQAChecker+;
#pragma link C++ class AliMUONTriggerQAChecker+;

#pragma link C++ class AliMUONTriggerUtilities+;
#endif


