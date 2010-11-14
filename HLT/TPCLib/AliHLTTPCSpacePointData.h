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

  static UInt_t GetSlice( UInt_t Id )  { return (Id>>25)&0x3F; }
  static UInt_t GetPatch( UInt_t Id )  { return (Id>>22)&0x7; }
  static UInt_t GetNumber( UInt_t Id ) { return Id&0x003FFFFF; }
  static UInt_t GetID( UInt_t Slice, UInt_t Patch, UInt_t Number ){
    return ((Slice&0x3F)<<25)+((Patch&0x7)<<22) + (Number&0x003FFFFF);
  }

  void SetID( UInt_t Slice, UInt_t Patch, UInt_t Number ){
    fID = GetID(Slice, Patch,Number);
  }
  UInt_t GetSlice() const { return GetSlice(fID); }
  UInt_t GetPatch() const { return GetPatch(fID); }
  UInt_t GetNumber() const { return GetNumber(fID); }

};
typedef struct AliHLTTPCSpacePointData AliHLTTPCSpacePointData;


#endif /* SPACEPOINTDATA_H */
