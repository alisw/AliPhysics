#ifndef ALIHLTITSTRACKDATA_H
#define ALIHLTITSTRACKDATA_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliHLTTPCRootTypes.h"
#include "AliExternalTrackParam.h"

/**
 * @struct AliHLTTPCTrackSegmentData
 * Primitive data exchange structure for TPC tracks.
 *
 * @ingroup alihlt_tpc_datastructs
 */
struct AliHLTITSTrackData
{
  AliExternalTrackParam fTrackParam; // track parameters in ITS
  Int_t fTPCId; // Id of the corresponding TPC track
  Int_t fClusterIds[6]; // Id's of the associated ITS clusters
};

typedef struct AliHLTITSTrackData AliHLTITSTrackData;

#endif /* _ALIHLTITSTRACKDATA_H_ */
