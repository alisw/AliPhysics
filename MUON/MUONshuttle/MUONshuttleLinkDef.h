/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file MUONshuttleLinkDef.h
/// \brief The CINT link definitions for \ref shuttle

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliMUONPreprocessor+;
#pragma link C++ class AliMUONTrackerPreprocessor+;
#pragma link C++ class AliMUONVSubprocessor+;
#pragma link C++ class AliMUONHVSubprocessor+;
#pragma link C++ class AliMUONLVSubprocessor+;
#pragma link C++ class AliMUONPedestalSubprocessor+;
#pragma link C++ class AliMUONGMSSubprocessor+;
#pragma link C++ class AliMUONTriggerSubprocessor+;
#pragma link C++ class AliMUONTriggerDCSSubprocessor+;
#pragma link C++ class AliMUONTriggerPreprocessor+;
#pragma link C++ class AliMUONOccupancySubprocessor+;
#pragma link C++ class AliMUONBusPatchEvolutionSubprocessor+;
#pragma link C++ class AliMUONConfigSubprocessor+;

#endif
