#ifndef ALIHLTPHOSCLUSTERDATASTRUCT
#define ALIHLTPHOSCLUSTERDATASTRUCT

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Øystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "Rtypes.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTDataTypes.h"

struct AliHLTPHOSClusterDataStruct
{
  Float_t fClusterEnergy;
  Float_t fLocalPositionPtr[2];
  AliHLTUInt8_t fPHOSModule;
};

#endif
