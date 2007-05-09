
#ifndef ALIHLTPHOSPHYSICSDEFINITIONS_H
#define ALIHLTPHOSPHYSICSDEFINITIONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTPHOSPhysicsDefinitions
 */

#include "AliHLTDataTypes.h"
#include "Rtypes.h"
#include "TROOT.h"
#include "TH1F.h"

class AliHLTPHOSPhysicsDefinitions
    {
    public:
      static const AliHLTComponentDataType fgkAliHLTClusterDataType;  
      static const AliHLTComponentDataType fgkAliHLTHistDataType;
      static const AliHLTComponentDataType fgkAliHLTSpectrumDataType;
    };

#endif

