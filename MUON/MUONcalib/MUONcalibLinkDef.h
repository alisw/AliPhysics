/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file MUONcalibLinkDef.h
/// \brief The CINT link definitions for \ref calib 


#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliMUONVCalibParam+;
#pragma link C++ class AliMUONCalibParamND+;
#pragma link C++ class AliMUONCalibParamNF+;
#pragma link C++ class AliMUONCalibParamNI+;
#pragma link C++ class AliMUONCalibrationData+;
#pragma link C++ class AliMUONTriggerEfficiencyCells+;
#pragma link C++ class AliMUONTriggerChamberEfficiency+;
#pragma link C++ class AliMUONTriggerLut+;
#pragma link C++ class AliMUONRegionalTriggerConfig+;
#pragma link C++ class AliMUONTriggerCrateConfig+;
#pragma link C++ class AliMUONGlobalCrateConfig+;
#pragma link C++ class AliMUONTriggerDisplay+;
#pragma link C++ class AliMUON2DStoreValidator+;
#pragma link C++ class AliMUONTriggerIO+;
#pragma link C++ class AliMUONTrackerIO+;
#pragma link C++ class AliMUONVTrackerData+;
#pragma link C++ class AliMUONSparseHisto+;
#pragma link C++ class AliMUONTrackerData-;
#pragma link C++ class AliMUONPedestal+;
#pragma link C++ class AliMUONErrorCounter+;
#pragma link C++ class AliMUONRejectList+;
#pragma link C++ class AliMUONTriggerScalers+;
#pragma link C++ class AliMUONBusPatchEvolution+;

#endif

