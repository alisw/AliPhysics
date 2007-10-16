// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTRDDEFINITIONS_H
#define ALIHLTTRDDEFINITIONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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
  
  static const AliHLTComponentDataType fgkDDLRawDataType; // Raw Data
  static const AliHLTComponentDataType fgkClusterDataType; // TRD Cluster Data
  static const AliHLTComponentDataType fgkTRDSATracksDataType; // Stand Alone tracks
  static const AliHLTComponentDataType fgkTRDSAEsdDataType; // Stand Alone tracks
  static const AliHLTComponentDataType fgkMCMtrackletDataType; // MCM tracklet Data
  static const AliHLTComponentDataType fgkMCMcalibrationDataType; // MCM Calibration data
  static const AliHLTComponentDataType fgkCalibrationDataType; // Calibration with TRDtracks

  ClassDef(AliHLTTRDDefinitions, 0)
    
};

#endif
