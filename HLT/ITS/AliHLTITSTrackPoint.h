#ifndef ALIHLTITSTRACKPOINT_H
#define ALIHLTITSTRACKPOINT_H

#include "AliHLTDataTypes.h"
#include "AliTrackPointArray.h"

// struct to hold the information on the space points
struct AliHLTITSTrackPoint 
{
  float fXYZ[3]; // coordinates
  float	fCov[6];	// Cov matrix
  float	fCharge;	// Cluster charge in arbitrary units
  float	fChargeRatio;	// Charge ratio in SSD
  int	fClusterType;	// Cluster Type (encoded info on size and shape)
  float	fDriftTime;	// Drift time in SDD (in ns)  
  unsigned short	fVolumeID;	// Volume ID
  void Reset();
  AliTrackPoint GetAliTrackPoint() const;
};
 
inline void AliHLTITSTrackPoint::Reset()
{
  //reset
  for( int i=0; i<3; i++) fXYZ[i] = 0;
  for( int i=0; i<6; i++) fCov[i] = 0;
  fCharge = 0;
  fChargeRatio = 0;
  fClusterType = 0;
  fVolumeID = 0;
}

inline AliTrackPoint AliHLTITSTrackPoint::GetAliTrackPoint() const
{
  return AliTrackPoint( fXYZ, fCov, fVolumeID, fCharge, fDriftTime, fChargeRatio, fClusterType);
}

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
