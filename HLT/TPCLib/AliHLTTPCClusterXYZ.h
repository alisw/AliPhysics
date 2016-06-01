// $Id$

#ifndef ALIHLTTPCCLUSTERXYZ_H
#define ALIHLTTPCCLUSTERXYZ_H

#include "Rtypes.h"

/**
 * @struct AliHLTTPCClusterXYZ
 * Primitive data exchange structure for XYZ (transformed) coordinates of TPC clusters.
 * Field fRawID points to corresponding AliHLTTPCRawCluster, where the original data is stored.
 *
 * @ingroup alihlt_tpc_datastructs 
 */

struct AliHLTTPCClusterXYZ{

  Float_t fX;            // X coordinate in local coordinates of HLT slice ( == inner + outer parts of TPC sector == 6 readout partitions )
  Float_t fY;            // Y coordinate in local coordinates
  Float_t fZ;            // Z coordinate in local coordinates
  UInt_t  fRawClusterID; // contains slice, partition, and index withing the partition of corresponding raw cluster

  void SetX( Float_t x )             { fX=x; }
  void SetY( Float_t y )             { fY=y; }
  void SetZ( Float_t z )             { fZ=z; }
  void SetRawClusterID( UInt_t id )  { fRawClusterID=id; }

  Float_t GetX()               const  { return fX; }
  Float_t GetY()               const  { return fY; }
  Float_t GetZ()               const  { return fZ; }
  UInt_t  GetRawClusterID()    const  { return fRawClusterID; }

  AliHLTTPCClusterXYZ() : fX(0.), fY(0.), fZ(0.), fRawClusterID(0) {}

};
typedef struct AliHLTTPCClusterXYZ AliHLTTPCClusterXYZ;


struct AliHLTTPCClusterXYZData
{
  UInt_t fCount;   // number of clusters
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTPCClusterXYZ  fClusters[1]; // array of clusters  
#else
  AliHLTTPCClusterXYZ  fClusters[0]; // array of clusters 
#endif
};
typedef struct AliHLTTPCClusterXYZData AliHLTTPCClusterXYZData;

#endif /* ALIHLTTPCCLUSTERXYZ_H */
