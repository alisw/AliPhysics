//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDDEFINITIONS_H
#define ALIHLTTRDDEFINITIONS_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/* AliHLTTRDDefinitions
 */

#include "AliHLTDataTypes.h"
#include "Rtypes.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  The HLT definitions for TRD                                              //  
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


class AliHLTTRDDefinitions
{

public:
  AliHLTTRDDefinitions();
  virtual ~AliHLTTRDDefinitions();
  
  static const AliHLTComponentDataType fgkClusterDataType; // TRD Cluster Data
  static const AliHLTComponentDataType fgkTRDSATracksDataType; // Stand Alone tracks
  static const AliHLTComponentDataType fgkTRDSAEsdDataType; // Stand Alone tracks
  static const AliHLTComponentDataType fgkMCMtrackletDataType; // MCM tracklet Data
  static const AliHLTComponentDataType fgkMCMcalibrationDataType; // MCM Calibration data
  static const AliHLTComponentDataType fgkCalibrationDataType; // Calibration with TRDtracks
  static const AliHLTComponentDataType fgkEORCalibrationDataType;//Calibration end of run 

  ClassDef(AliHLTTRDDefinitions, 0)
    
};

#endif

