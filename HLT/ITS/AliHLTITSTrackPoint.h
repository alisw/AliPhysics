#ifndef ALIHLTITSTRACKPOINT_H
#define ALIHLTITSTRACKPOINT_H

#include "AliHLTDataTypes.h"

// struct to hold the information on the space points
struct AliHLTITSTrackPoint {
  float fX[3];
  unsigned short fVolumeId; 
};

struct AliHLTITSTrackPointData {
  AliHLTUInt32_t fCount; // number of space points
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTITSTrackPoint fPoints[1]; // array of space points
#else
  AliHLTITSTrackPoint fPoints[0]; // array of space points
#endif
};

typedef struct AliHLTITSTrackPointData AliHLTITSTrackPointData;


#endif
