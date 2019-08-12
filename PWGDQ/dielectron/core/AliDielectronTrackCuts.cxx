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

///////////////////////////////////////////////////////////////////////////
//                Dielectron TrackCuts                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include <TMath.h>

#include "AliDielectronTrackCuts.h"
#include "AliDielectronClusterCuts.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"

ClassImp(AliDielectronTrackCuts)

AliDielectronTrackCuts::AliDielectronTrackCuts() :
  AliAnalysisCuts(),
  fV0DaughterCut(0),
  fNegateV0DauterCut(kFALSE),
  fITSclusterBitMap(0),
  fITSclusterCutType(kOneOf),
  fSelectGlobalTrack(kFALSE),
  fRequireITSRefit(kFALSE),
  fRequireTPCRefit(kFALSE),
  fTPCNclRobustCut(-1),
  fTPCcrossedOverFindable(-1.),
  fAODFilterBit(kSwitchOff),
  fWaiveITSNcls(-1),
  fRequireTRDUpdate(kFALSE),
  fRequireCaloClusterMatch(kFALSE),
  fClusterMatchCaloType(AliDielectronClusterCuts::kAny)
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
  fITSclusterBitMap(0),
  fITSclusterCutType(kOneOf),
  fSelectGlobalTrack(kFALSE),
  fRequireITSRefit(kFALSE),
  fRequireTPCRefit(kFALSE),
  fTPCNclRobustCut(-1),
  fTPCcrossedOverFindable(-1.),
  fAODFilterBit(kSwitchOff),
  fWaiveITSNcls(-1),
  fRequireTRDUpdate(kFALSE),
  fRequireCaloClusterMatch(kFALSE),
  fClusterMatchCaloType(AliDielectronClusterCuts::kAny)
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

  // use filter bit to speed up the AOD analysis (track pre-filter)
  // relevant filter bits are:
  // kTPCqual==1             -> TPC quality cuts
  // kTPCqualSPDany==4       -> + SPD any
  // kTPCqualSPDanyPIDele==8 -> + nSigmaTPCele +-3 (inclusion)

  if (track->IsA()==AliAODTrack::Class()) {
    if(fSelectGlobalTrack) if(((AliAODTrack*)track)->GetID() < 0) accept=kFALSE;
    if(fAODFilterBit!=kSwitchOff) accept*=((AliAODTrack*)track)->TestFilterBit(fAODFilterBit);
  }

  if (fV0DaughterCut) {
    Bool_t isV0=track->TestBit(BIT(fV0DaughterCut));
    if (fNegateV0DauterCut) isV0=!isV0;
    accept*=isV0;
  }

  //ESD track cut like ITS cluster cut
  for (Int_t i=0;i<3;++i){
    Bool_t layer1=TESTBIT(vtrack->GetITSClusterMap(),i*2);
    Bool_t layer2=TESTBIT(vtrack->GetITSClusterMap(),i*2+1);
    accept*=CheckITSClusterRequirement(fCutClusterRequirementITS[i], layer1, layer2);
  }

  //more flexible ITS cluster cut
  if (fITSclusterBitMap) accept*=CheckITSClusterCut(vtrack->GetITSClusterMap());

  //different its cluster cut
  if (fWaiveITSNcls > -1) {
    Int_t nITScls      = 0;
    Int_t requiredNcls = 7;
    for(Int_t i=5; i>=0; i--) {
      if(TESTBIT(vtrack->GetITSClusterMap(),i)) {
	nITScls++;
	requiredNcls=6-fWaiveITSNcls-i;
      }
    }
    accept*=(requiredNcls<=nITScls);
  }

  //its and tpc refit
  if (fRequireITSRefit) accept*=(vtrack->GetStatus()&AliVTrack::kITSrefit)>0;
  if (fRequireTPCRefit) accept*=(vtrack->GetStatus()&AliVTrack::kTPCrefit)>0;

  Int_t nclr=0;
  if (fTPCNclRobustCut>0){
    nclr=TMath::Nint(vtrack->GetTPCClusterInfo(2,1));
    accept*=(nclr>fTPCNclRobustCut);
  }
  if (fTPCcrossedOverFindable > 0.) {
    if(fTPCNclRobustCut<=0) nclr=TMath::Nint(vtrack->GetTPCClusterInfo(2,1));
    Int_t tpcNclsF = vtrack->GetTPCNclsF();
    accept*=(tpcNclsF); //ESDtrackCut would return here true
    if (tpcNclsF != 0) {//'accept' already negated above in this case above
      accept*=(((Double_t)nclr/(Double_t)vtrack->GetTPCNclsF()) >= fTPCcrossedOverFindable);
    }
  }

  // TRD update
  if (fRequireTRDUpdate) accept*=(vtrack->GetStatus()&AliVTrack::kTRDupdate)>0;

  // calo cluster-track match
  if (fRequireCaloClusterMatch) {
    Int_t fCaloIndex = vtrack->GetEMCALcluster();
    if (fCaloIndex!=AliVTrack::kEMCALNoMatch) {
      if (fClusterMatchCaloType==AliDielectronClusterCuts::kEMCal)  accept*=vtrack->IsEMCAL();
      if (fClusterMatchCaloType==AliDielectronClusterCuts::kPHOS)   accept*=vtrack->IsPHOS();
      if (fClusterMatchCaloType==AliDielectronClusterCuts::kAny)    accept*=kTRUE;
    } else {
      accept*=kFALSE;
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

//______________________________________________
Bool_t AliDielectronTrackCuts::CheckITSClusterCut(UChar_t itsBits) const
{
  // check the its cluster cut
  switch (fITSclusterCutType){
  case kOneOf:   return itsBits & fITSclusterBitMap;
  case kAtLeast: return (itsBits & fITSclusterBitMap)==fITSclusterBitMap;
  case kExact:   return (itsBits==fITSclusterBitMap);
  }
  return kTRUE;
}
