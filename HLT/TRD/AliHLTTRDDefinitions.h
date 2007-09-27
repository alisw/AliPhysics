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

      static const AliHLTComponentDataType fgkDDLRawDataType;
      static const AliHLTComponentDataType fgkClusterDataType;
      static const AliHLTComponentDataType fgkTRDSATracksDataType; // Stand Alone tracks
      static const AliHLTComponentDataType fgkTRDSAEsdDataType; // Stand Alone tracks
      static const AliHLTComponentDataType fgkMCMtrackletDataType;
      static const AliHLTComponentDataType fgkMCMcalibrationDataType;
      
      ClassDef(AliHLTTRDDefinitions, 0)

};

#endif
