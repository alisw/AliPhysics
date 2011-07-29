/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file MUONcalignLinkDef.h
/// \brief The CINT link definitions for \ref calign 


#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliMUONClusterInfo+;
#pragma link C++ class AliMUONPadInfo+;
#pragma link C++ class AliMUONAlignment+;
#pragma link C++ class AliMUONAlignmentClusterRecord+;
#pragma link C++ class AliMUONAlignmentTrackRecord+;
#pragma link C++ class AliMUONAlignmentTask+;
#pragma link C++ class AliMUONReAlignTask+;
#pragma link C++ class AliMUONChamberCalibrationTask+;

#endif

