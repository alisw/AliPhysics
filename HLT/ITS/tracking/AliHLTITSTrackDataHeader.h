#ifndef _ALIHLTITSTRACKDATAHEADER_HPP_
#define _ALIHLTITSTRACKDATAHEADER_HPP_

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliHLTDataTypes.h"
#include "AliHLTITSTrackData.h"

/**
 * @struct AliHLTITSTrackDataHeader
 * Primitive data exchange structure for ITS tracks.
 *
 * @ingroup alihlt_tpc_datastructs
 */
struct AliHLTITSTrackDataHeader
{
  AliHLTUInt32_t fTrackletCnt;
#ifndef __SUNPRO_CC
  AliHLTITSTrackData fTracks[];
#else
  AliHLTITSTrackData fTracks[1];
#endif
};

typedef struct AliHLTITSTrackDataHeader AliHLTITSTrackDataHeader;

#endif // _ALIHLTITSTRACKDATAHEADER_HPP_
