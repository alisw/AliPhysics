// $Id: AliHLTExternalTrackParam.h 33222 2009-06-26 23:45:52Z sgorbuno $
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef _ALIHLTTRACKMCDATA_H_
#define _ALIHLTTRACKMCDATA_H_

#include "AliHLTDataTypes.h"

/**
 * @struct AliHLTTrackMCData
 * This in a struct for MC track labels
 * @ingroup alihlt_component_datatypes
 */
struct AliHLTTrackMCLabel
{
  Int_t fTrackID;
  Int_t fMCLabel;
};

typedef struct AliHLTTrackMCLabel AliHLTTrackMCLabel;

struct AliHLTTrackMCData {
  AliHLTUInt32_t fCount;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTrackMCLabel fLabels[1];
#else
  AliHLTTrackMCLabel fLabels[];
#endif
};

typedef struct AliHLTTrackMCData AliHLTTrackMCData;

#endif
