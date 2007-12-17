#include <Rtypes.h>

struct AliHLTTPCSpacePointData
{
  //#ifdef do_mc
  //Int_t fTrackID[3];
  //#endif
  Float_t fX;  //==fPadRow in local system
  Float_t fY;  
  Float_t fZ;
  UInt_t fID;  //contains slice patch and number
  UChar_t fPadRow;
  Float_t fSigmaY2; //error (former width) of the clusters
  Float_t fSigmaZ2; //error (former width) of the clusters
  UInt_t fCharge;
  Bool_t fUsed;     // only used in AliHLTTPCDisplay 
  Int_t fTrackN;    // only used in AliHLTTPCDisplay 
};

struct AliHLTTPCClusterData
{
  Int_t fSpacePointCnt;
  AliHLTTPCSpacePointData fSpacePoints[1];
};

void hlt_structs()
{}
