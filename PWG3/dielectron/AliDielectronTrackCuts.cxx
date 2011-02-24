/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id$ */


///////////////////////////////////////////////////////////////////////////
//                Dielectron TrackCuts                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////



#include "AliDielectronTrackCuts.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"

ClassImp(AliDielectronTrackCuts)

AliDielectronTrackCuts::AliDielectronTrackCuts() :
  AliAnalysisCuts(),
  fV0DaughterCut(0),
  fNegateV0DauterCut(kFALSE),
  fRequireITSRefit(kFALSE),
  fRequireTPCRefit(kFALSE),
  fTPCNclRobustCut(-1)
{
  //
  // Default Constructor
  //

  for (Int_t i = 0; i < 3; i++)
    fCutClusterRequirementITS[i] = kOff;
  
}

//______________________________________________
AliDielectronTrackCuts::AliDielectronTrackCuts(const char* name, const char* title) :
  AliAnalysisCuts(name, title),
  fV0DaughterCut(0),
  fNegateV0DauterCut(kFALSE),
  fRequireITSRefit(kFALSE),
  fRequireTPCRefit(kFALSE),
  fTPCNclRobustCut(-1)
{
  //
  // Named Constructor
  //

  for (Int_t i = 0; i < 3; i++)
    fCutClusterRequirementITS[i] = kOff;
  
}

//______________________________________________
AliDielectronTrackCuts::~AliDielectronTrackCuts()
{
  //
  // Default Destructor
  //
  
}

//______________________________________________
Bool_t AliDielectronTrackCuts::IsSelected(TObject* track)
{
  //
  // Apply configured cuts
  //

  AliVTrack *vtrack=dynamic_cast<AliVTrack*>(track);
  if (!vtrack) return kFALSE;
  
  Bool_t accept=kTRUE;
  if (fV0DaughterCut) {
    Bool_t isV0=track->TestBit(BIT(fV0DaughterCut));
    if (fNegateV0DauterCut) isV0=!isV0;
    accept*=isV0;
  }

  for (Int_t i=0;i<3;++i){
    Bool_t layer1=TESTBIT(vtrack->GetITSClusterMap(),i*2);
    Bool_t layer2=TESTBIT(vtrack->GetITSClusterMap(),i*2+1);
    accept*=CheckITSClusterRequirement(fCutClusterRequirementITS[i], layer1, layer2);
  }

  if (fRequireITSRefit) accept*=(vtrack->GetStatus()&kITSrefit)>0;
  if (fRequireTPCRefit) accept*=(vtrack->GetStatus()&kTPCrefit)>0;

  if (fTPCNclRobustCut>0){
    AliESDtrack *tr=dynamic_cast<AliESDtrack*>(track);
    if (tr){
      Int_t nclr=TMath::Nint(tr->GetTPCClusterInfo(2,1));
      accept*=(nclr>fTPCNclRobustCut);
    }
  }
  return accept;
}

//______________________________________________
void AliDielectronTrackCuts::SetV0DaughterCut(AliPID::EParticleType type, Bool_t negate/*=kFALSE*/)
{
  //
  // Set V0 Daughter cut bit
  //
  const Int_t bitMap[5] = {14, -1, 15, -1, 16}; // convert the AliPID to bit positions
  fV0DaughterCut=bitMap[type];
  fNegateV0DauterCut=negate;
}

//____________________________________________________________________
Bool_t AliDielectronTrackCuts::CheckITSClusterRequirement(ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2) const
{
  // checks if the cluster requirement is fullfilled (in this case: return kTRUE)
  
  switch (req)
  {
  case kOff:        return kTRUE;
  case kNone:       return !clusterL1 && !clusterL2;
  case kAny:        return clusterL1 || clusterL2;
  case kFirst:      return clusterL1;
  case kOnlyFirst:  return clusterL1 && !clusterL2;
  case kSecond:     return clusterL2;
  case kOnlySecond: return clusterL2 && !clusterL1;
  case kBoth:       return clusterL1 && clusterL2;
  }
  
  return kFALSE;
}

