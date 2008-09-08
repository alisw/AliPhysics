// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <Rtypes.h>
#include "TEveElement.h"
#include "TPolyMarker3D.h"

struct AliHLTTPCSpacePointData
{
  Float_t fX;       // == fPadRow in local system
  Float_t fY;
  Float_t fZ;
  UInt_t fID;       //contains slice patch and number
  UChar_t fPadRow;
  Float_t fSigmaY2; //error (former width) of the clusters
  Float_t fSigmaZ2; //error (former width) of the clusters
  UInt_t fCharge;
  UInt_t fMaxQ;     // QMax of cluster
  Bool_t fUsed;     // only used in AliHLTTPCDisplay
  Int_t fTrackN;    // only used in AliHLTTPCDisplay
};

struct AliHLTTPCClusterData
{
  Int_t fSpacePointCnt;
  AliHLTTPCSpacePointData fSpacePoints[1];
};


struct AliHLTMUONTrackerRawData{
  UInt_t fDataId;
  UShort_t fADC;
};

struct AliHLTMUONTrackerMappingData{
  Int_t fDetElemId;
  Int_t fIX,fIY;
  Int_t fCath;
  Float_t fRealX,fRealY,fRealZ;
  Int_t fADC;
  Float_t fPed;
  Float_t fSigma;
};

struct AliHLTMUONTriggerPointData{
  Int_t fDetElemId;
  Float_t fX,fY,fZ ;
};

struct AliHLTMUONTriggerMappingData{
  AliHLTMUONTriggerPointData fLut[8][16][4][2][16];
};

void hlt_structs()
{}
