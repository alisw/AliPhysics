// @(#) $Id$

#ifndef SPACEPOINTDATA_H
#define SPACEPOINTDATA_H

#include "AliHLTRootTypes.h"
struct AliHLTSpacePointData{
#ifdef do_mc
  Int_t fTrackID[3];
#endif
  Float_t fX;  //==fPadRow in local system
  Float_t fY;  
  Float_t fZ;
  UInt_t fID;  //contains slice patch and number
  UChar_t fPadRow;
  Float_t fSigmaY2; //error (former width) of the clusters
  Float_t fSigmaZ2; //error (former width) of the clusters
  UInt_t fCharge;
};
typedef struct AliHLTSpacePointData AliHLTSpacePointData;
typedef AliHLTSpacePointData AliL3SpacePointData;


#endif /* SPACEPOINTDATA_H */
