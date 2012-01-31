/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// AliStarTrackCuts:
// A track cut class for the AliStarTrack
//
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include <limits.h>
#include <float.h>
#include "TNamed.h"
#include "AliStarTrack.h"
#include "AliStarTrackCuts.h"

ClassImp(AliStarTrackCuts)

//-----------------------------------------------------------------------
AliStarTrackCuts::AliStarTrackCuts():
  TNamed(),
  fCutID(kFALSE),
  fIDMax(INT_MAX),
  fIDMin(INT_MIN),
  fCutCharge(kFALSE),
  fChargeMax(INT_MAX),
  fChargeMin(INT_MIN),
  fCutEta(kFALSE),
  fEtaMax(FLT_MAX),
  fEtaMin(-FLT_MAX),
  fCutPhi(kFALSE),
  fPhiMax(FLT_MAX),
  fPhiMin(-FLT_MAX),
  fCutPt(kFALSE),
  fPtMax(FLT_MAX),
  fPtMin(-FLT_MAX),
  fCutDCA(kFALSE),
  fDCAMax(FLT_MAX),
  fDCAMin(-FLT_MAX),
  fCutNHits(kFALSE),
  fNHitsMax(INT_MAX),
  fNHitsMin(INT_MIN),
  fCutNHitsFit(kFALSE),
  fNHitsFitMax(INT_MAX),
  fNHitsFitMin(INT_MIN),
  fCutNHitsPoss(kFALSE),
  fNHitsPossMax(INT_MAX),
  fNHitsPossMin(INT_MIN),
  fCutNHitsDedx(kFALSE),
  fNHitsDedxMax(INT_MAX),
  fNHitsDedxMin(INT_MIN),
  fCutdEdx(kFALSE),
  fdEdxMax(FLT_MAX),
  fdEdxMin(-FLT_MAX),
  fCutNSigElect(kFALSE),
  fNSigElectMax(FLT_MAX),
  fNSigElectMin(-FLT_MAX),
  fCutNSigPi(kFALSE),
  fNSigPiMax(FLT_MAX),
  fNSigPiMin(-FLT_MAX),
  fCutNSigK(kFALSE),
  fNSigKMax(FLT_MAX),
  fNSigKMin(-FLT_MAX),
  fCutNSigProton(kFALSE),
  fNSigProtonMax(FLT_MAX),
  fNSigProtonMin(-FLT_MAX),
  fCutFitRatio(kFALSE),
  fFitRatioMax(FLT_MAX),
  fFitRatioMin(FLT_MIN)
{
  //ctor 
}

////-----------------------------------------------------------------------
//AliStarTrackCuts::AliStarTrackCuts(const AliStarTrackCuts& that):
//  TNamed(),
//  fCutID(that.fCutID),
//  fIDMax(that.fIDMax),
//  fIDMin(that.fIDMin),
//  fCutCharge(that.fCutCharge),
//  fChargeMax(that.fChargeMax),
//  fChargeMin(that.fChargeMin),
//  fCutEta(that.fCutEta),
//  fEtaMax(that.fEtaMax),
//  fEtaMin(that.fEtaMin),
//  fCutPhi(that.fCutPhi),
//  fPhiMax(that.fPhiMax),
//  fPhiMin(that.fPhiMin),
//  fCutPt(that.fCutPt),
//  fPtMax(that.fPtMax),
//  fPtMin(that.fPtMin),
//  fCutDCA(that.fCutDCA),
//  fDCAMax(that.fDCAMax),
//  fDCAMin(that.fDCAMin),
//  fCutNHits(that.fCutNHits),
//  fNHitsMax(that.fNHitsMax),
//  fNHitsMin(that.fNHitsMin),
//  fCutNHitsFit(that.fCutNHitsFit),
//  fNHitsFitMax(that.fNHitsFitMax),
//  fNHitsFitMin(that.fNHitsFitMin),
//  fCutNHitsPoss(that.fCutNHitsPoss),
//  fNHitsPossMax(that.fNHitsPossMax),
//  fNHitsPossMin(that.fNHitsPossMin),
//  fCutNHitsDedx(that.fCutNHitsDedx),
//  fNHitsDedxMax(that.fNHitsDedxMax),
//  fNHitsDedxMin(that.fNHitsDedxMin),
//  fCutdEdx(that.fCutdEdx),
//  fdEdxMax(that.fdEdxMax),
//  fdEdxMin(that.fdEdxMin),
//  fCutNSigElect(that.fCutNSigElect),
//  fNSigElectMax(that.fNSigElectMax),
//  fNSigElectMin(that.fNSigElectMin),
//  fCutNSigPi(that.fCutNSigPi),
//  fNSigPiMax(that.fNSigPiMax),
//  fNSigPiMin(that.fNSigPiMin),
//  fCutNSigK(that.fCutNSigK),
//  fNSigKMax(that.fNSigKMax),
//  fNSigKMin(that.fNSigKMin),
//  fCutNSigProton(that.fCutNSigProton),
//  fNSigProtonMax(that.fNSigProtonMax),
//  fNSigProtonMin(that.fNSigProtonMin),
//  fCutFitRatio(that.fCutFitRatio),
//  fFitRatioMax(that.fFitRatioMax),
//  fFitRatioMin(that.fFitRatioMin)
//{
//  //copy ctor 
//}
//
////-----------------------------------------------------------------------
//AliStarTrackCuts& AliStarTrackCuts::operator=(const AliStarTrackCuts& that)
//{
//  //assignment
//  fCutID=that.fCutID;
//  fIDMax=that.fIDMax;
//  fIDMin=that.fIDMin;
//  fCutCharge=that.fCutCharge;
//  fChargeMax=that.fChargeMax;
//  fChargeMin=that.fChargeMin;
//  fCutEta=that.fCutEta;
//  fEtaMax=that.fEtaMax;
//  fEtaMin=that.fEtaMin;
//  fCutPhi=that.fCutPhi;
//  fPhiMax=that.fPhiMax;
//  fPhiMin=that.fPhiMin;
//  fCutPt=that.fCutPt;
//  fPtMax=that.fPtMax;
//  fPtMin=that.fPtMin;
//  fCutDCA=that.fCutDCA;
//  fDCAMax=that.fDCAMax;
//  fDCAMin=that.fDCAMin;
//  fCutNHits=that.fCutNHits;
//  fNHitsMax=that.fNHitsMax;
//  fNHitsMin=that.fNHitsMin;
//  fCutNHitsFit=that.fCutNHitsFit;
//  fNHitsFitMax=that.fNHitsFitMax;
//  fNHitsFitMin=that.fNHitsFitMin;
//  fCutNHitsPoss=that.fCutNHitsPoss;
//  fNHitsPossMax=that.fNHitsPossMax;
//  fNHitsPossMin=that.fNHitsPossMin;
//  fCutNHitsDedx=that.fCutNHitsDedx;
//  fNHitsDedxMax=that.fNHitsDedxMax;
//  fNHitsDedxMin=that.fNHitsDedxMin;
//  fCutdEdx=that.fCutdEdx;
//  fdEdxMax=that.fdEdxMax;
//  fdEdxMin=that.fdEdxMin;
//  fCutNSigElect=that.fCutNSigElect;
//  fNSigElectMax=that.fNSigElectMax;
//  fNSigElectMin=that.fNSigElectMin;
//  fCutNSigPi=that.fCutNSigPi;
//  fNSigPiMax=that.fNSigPiMax;
//  fNSigPiMin=that.fNSigPiMin;
//  fCutNSigK=that.fCutNSigK;
//  fNSigKMax=that.fNSigKMax;
//  fNSigKMin=that.fNSigKMin;
//  fCutNSigProton=that.fCutNSigProton;
//  fNSigProtonMax=that.fNSigProtonMax;
//  fNSigProtonMin=that.fNSigProtonMin;
//  fCutFitRatio=that.fCutFitRatio;
//  fFitRatioMax=that.fFitRatioMax;
//  fFitRatioMin=that.fFitRatioMin;
//  
//  return *this;
//}

//----------------------------------------------------------------------- 
Bool_t AliStarTrackCuts::PassesCuts(const AliStarTrack *track) const
{
  //check is track passes cuts
  if(fCutID) {if (track->GetID() < fIDMin || track->GetID() > fIDMax ) return kFALSE;} //integer values: non inclusive bound!
  if(fCutCharge) {if (track->GetCharge() < fChargeMin || track->GetCharge() > fChargeMax ) return kFALSE;}
  if(fCutEta) {if (track->GetEta() < fEtaMin || track->GetEta() >= fEtaMax ) return kFALSE;}
  if(fCutPhi) {if (track->GetPhi() < fPhiMin || track->GetPhi() >= fPhiMax ) return kFALSE;}
  if(fCutPt) {if (track->GetPt() < fPtMin || track->GetPt() >= fPtMax ) return kFALSE;}
  if(fCutDCA) {if (track->GetDCA() < fDCAMin || track->GetDCA() >= fDCAMax ) return kFALSE;}
  if(fCutNHits) {if (track->GetNHits() < fNHitsMin || track->GetNHits() > fNHitsMax ) return kFALSE;}
  if(fCutNHitsFit) {if (track->GetNHitsFit() < fNHitsFitMin || track->GetNHitsFit() > fNHitsFitMax ) return kFALSE;}
  if(fCutNHitsPoss) {if (track->GetNHitsPoss() < fNHitsPossMin || track->GetNHitsPoss() > fNHitsPossMax ) return kFALSE;}
  if(fCutNHitsDedx) {if (track->GetNHitsDedx() < fNHitsDedxMin || track->GetNHitsDedx() > fNHitsDedxMax ) return kFALSE;}
  if(fCutdEdx) {if (track->GetdEdx() < fdEdxMin || track->GetdEdx() >= fdEdxMax ) return kFALSE;}
  if(fCutNSigElect) {if (track->GetNSigElect() < fNSigElectMin || track->GetNSigElect() > fNSigElectMax ) return kFALSE;}
  if(fCutNSigPi) {if (track->GetNSigPi() < fNSigPiMin || track->GetNSigPi() > fNSigPiMax ) return kFALSE;}
  if(fCutNSigK) {if (track->GetNSigK() < fNSigKMin || track->GetNSigK() > fNSigKMax ) return kFALSE;}
  if(fCutNSigProton) {if (track->GetNSigProton() < fNSigProtonMin || track->GetNSigProton() > fNSigProtonMax ) return kFALSE;}
  if(fCutFitRatio)
  {
    Int_t nhitsposs =  track->GetNHitsPoss();
    if (nhitsposs==0) return kFALSE;
    Float_t ratio = Float_t(track->GetNHitsFit()) / nhitsposs;
    if ( ratio < fFitRatioMin && ratio >= fFitRatioMax ) return kFALSE;
  }

  return kTRUE;
}

//----------------------------------------------------------------------- 
AliStarTrackCuts* AliStarTrackCuts::StandardCuts()
{
  //return a set of standard cuts, caller becomes owner
  AliStarTrackCuts* cuts = new AliStarTrackCuts();
  cuts->SetDCAMin(0.0);     // cm
  cuts->SetDCAMax(3.0);      
  cuts->SetPtMin(0.15);      // GeV
  cuts->SetPtMax(8.0);
  cuts->SetEtaMin(-1.1);
  cuts->SetEtaMax(1.1);
  cuts->SetFitRatioMin(0.52);  // Number of hits over number of hits possible
  cuts->SetFitRatioMax(2.0);      //should not exceed 1.0 but who knows...
  cuts->SetNHitsPossMin(5);    // Don't bother to fit tracks if # possible hits is too low
  cuts->SetNHitsPossMax(100);
  cuts->SetNHitsFitMin(15);    // 15 is typical but sometimes goes as high as 25
  cuts->SetNHitsFitMax(100);   // 45 pad rows in the TPC and so anything bigger than 45+Silicon is infinite
  return cuts;
}
