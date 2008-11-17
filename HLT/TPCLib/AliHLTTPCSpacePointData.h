// @(#) $Id$
// Original: AliHLTSpacePointData.h,v 1.4 2003/07/27 21:02:09 loizides 

#ifndef SPACEPOINTDATA_H
#define SPACEPOINTDATA_H

#include "AliHLTTPCRootTypes.h"

/**
 * @struct AliHLTTPCSpacePointData
 * Primitive data exchange structure for TPC clusters.
 * Together with the AliHLTTPCClusterDataFormat this defines
 * the output of the TPC online Cluster Finder.
 *
 * To translate between local coordinates, global coordinates, and row-pad-time coordinates 
 * one cann use AliHLTTPCTransform. 
 *
 * See AliHLTTPCClusterFinder::WriteClusters() for example. 
 *
 * @ingroup alihlt_tpc_datastructs
 */
struct AliHLTTPCSpacePointData{
  Float_t fX;       // X coordinate in local coordinates
  Float_t fY;       // Y coordinate in local coordinates
  Float_t fZ;       // Z coordinate in local coordinates
  UInt_t fID;       // contains slice patch and number
  UChar_t fPadRow;  // Pad row number
  Float_t fSigmaY2; // error (former width) of the clusters
  Float_t fSigmaZ2; // error (former width) of the clusters
  UInt_t fCharge;   // total charge of cluster
  UInt_t fQMax;     // QMax of cluster
  Bool_t fUsed;     // only used in AliHLTTPCDisplay 
  Int_t fTrackN;    // only used in AliHLTTPCDisplay 
};
typedef struct AliHLTTPCSpacePointData AliHLTTPCSpacePointData;


#endif /* SPACEPOINTDATA_H */
