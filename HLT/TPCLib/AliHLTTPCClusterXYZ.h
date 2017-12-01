// $Id$

#ifndef ALIHLTTPCCLUSTERXYZ_H
#define ALIHLTTPCCLUSTERXYZ_H

#include "Rtypes.h"

/**
 * @struct AliHLTTPCClusterXYZ
 *
 * Primitive data exchange structure for XYZ (transformed) coordinates of TPC clusters.
 * There is one-to-one correspondance between XYZ arrays and raw cluster arrays for each TPC partition, no extra indexing is needed.
 *
 * In the case i-th raw cluster is currupted and can not be transformed, its XYZ coordinates will be nevertheless present in the array at position i, 
 * but set to {0,0,0}
 *
 * @ingroup alihlt_tpc_datastructs 
 */

struct AliHLTTPCClusterXYZ{

  Float_t fX;            // X coordinate in local coordinates of HLT slice. ( one slice == geometrical TPC sector == inner + outer offline sectors == 6 readout partitions )
  Float_t fY;            // Y coordinate in local coordinates of HLT slice.
  Float_t fZ;            // Z coordinate in local coordinates of HLT slice.

  void Set( Float_t x, Float_t y, Float_t z ){ fX=x; fY=y; fZ=z; }

  void SetX( Float_t x )  { fX=x; }
  void SetY( Float_t y )  { fY=y; }
  void SetZ( Float_t z )  { fZ=z; }

  Float_t GetX()  const  { return fX; }
  Float_t GetY()  const  { return fY; }
  Float_t GetZ()  const  { return fZ; }

  AliHLTTPCClusterXYZ() : fX(0.), fY(0.), fZ(0.) {}

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
