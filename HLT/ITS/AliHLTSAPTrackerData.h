// $Id$
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef ALIHLTITSSAPTRACKERDATA_H
#define ALIHLTITSSAPTRACKERDATA_H

#include "AliHLTDataTypes.h"
#include "AliFlatExternalTrackParam.h"
#include "AliHLTStdIncludes.h"

struct AliHLTITSSAPTrackerData
{ 
  AliFlatExternalTrackParam paramOut;
  AliFlatExternalTrackParam paramInw;
  float chi2;
  short ncl;
  int   label;
};

typedef struct AliHLTITSSAPTrackerData AliHLTITSSAPTrackerData;

struct AliHLTITSSAPTrackerDataContainer {
  AliHLTUInt32_t fCount; // number of tracklets
  AliHLTUInt32_t fNSPDtracklets;
  AliHLTUInt32_t fNclusters[6];
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTITSSAPTrackerData fTracks[1]; // array of tracklets
#else
  AliHLTITSSAPTrackerData fTracks[0]; // array of tracklets
#endif
};

typedef struct AliHLTITSSAPTrackerDataContainer AliHLTITSSAPTrackerDataContainer;

#endif
