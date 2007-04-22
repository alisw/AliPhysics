// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTRDDEFINITIONS_H
#define ALIHLTTRDDEFINITIONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTRDDefinitions
 */

#include "AliHLTDataTypes.h"
#include "TObject.h"

class AliHLTTRDDefinitions
    {
    public:

      static const AliHLTComponentDataType gkDDLRawDataType;
      static const AliHLTComponentDataType gkClusterDataType;
      static const AliHLTComponentDataType gkTRDSATracksDataType; // Stand Alone tracks
      static const AliHLTComponentDataType gkMCMtrackletDataType;
      static const AliHLTComponentDataType gkMCMcalibrationDataType;
      
      ClassDef(AliHLTTRDDefinitions, 0)

    };

typedef struct AliTRDDummyData
{
  char str[50];
  long int val;
};

#endif
