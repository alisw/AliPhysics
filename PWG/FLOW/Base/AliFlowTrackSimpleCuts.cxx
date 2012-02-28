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

// AliFlowTrackSimpleCuts:
// A simple track cut class to the the AliFlowTrackSimple 
// for basic kinematic cuts
//
// author: N. van der Kolk (kolk@nikhef.nl)
// mods: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include <limits.h>
#include <float.h>
#include "TNamed.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "AliFlowTrackSimpleCuts.h"
#include "AliFlowTrackSimple.h"

ClassImp(AliFlowTrackSimpleCuts)

//-----------------------------------------------------------------------
AliFlowTrackSimpleCuts::AliFlowTrackSimpleCuts(const char* name):
  TNamed(name,name),
  fCutPt(kFALSE),
  fPtMax(FLT_MAX),
  fPtMin(-FLT_MAX),
  fCutEta(kFALSE),
  fEtaMax(FLT_MAX),
  fEtaMin(-FLT_MAX),
  fCutPhi(kFALSE),
  fPhiMax(FLT_MAX),
  fPhiMin(-FLT_MAX),
  fCutPID(kFALSE),
  fPID(0),
  fCutCharge(kFALSE),
  fCharge(0),
  fCutMass(kFALSE),
  fMassMax(FLT_MAX),
  fMassMin(-FLT_MAX)
{
  //constructor 
}

////-----------------------------------------------------------------------
//AliFlowTrackSimpleCuts::AliFlowTrackSimpleCuts(const AliFlowTrackSimpleCuts& someCuts):
//  TNamed(),
//  fCutPt(someCuts.fCutPt),
//  fPtMax(someCuts.fPtMax),
//  fPtMin(someCuts.fPtMin),
//  fCutEta(someCuts.fCutEta),
//  fEtaMax(someCuts.fEtaMax),
//  fEtaMin(someCuts.fEtaMin),
//  fCutPhi(someCuts.fCutPhi),
//  fPhiMax(someCuts.fPhiMax),
//  fPhiMin(someCuts.fPhiMin),
//  fCutPID(someCuts.fCutPID),
//  fPID(someCuts.fPID),
//  fCutCharge(someCuts.fCutCharge),
//  fCharge(someCuts.fCharge)
//{
//  //copy constructor 
//}
//
////-----------------------------------------------------------------------
//AliFlowTrackSimpleCuts& AliFlowTrackSimpleCuts::operator=(const AliFlowTrackSimpleCuts& someCuts)
//{
//  TNamed::operator=(someCuts);
//  fCutPt  = someCuts.fCutPt;
//  fPtMax  = someCuts.fPtMax;
//  fPtMin  = someCuts.fPtMin;
//  fCutEta = someCuts.fCutEta;
//  fEtaMax = someCuts.fEtaMax;
//  fEtaMin = someCuts.fEtaMin;
//  fCutPhi = someCuts.fCutPhi;
//  fPhiMax = someCuts.fPhiMax;
//  fPhiMin = someCuts.fPhiMin;
//  fCutPID = someCuts.fCutPID;
//  fPID    = someCuts.fPID;
//  fCutCharge = someCuts.fCutCharge;
//  fCharge = someCuts.fCharge;
//
//  return *this;
//}

//----------------------------------------------------------------------- 
Bool_t AliFlowTrackSimpleCuts::IsSelected(TObject* obj, Int_t)
{
  //check cuts
  TParticle* p = dynamic_cast<TParticle*>(obj);
  if (p) return PassesCuts(p);
  AliFlowTrackSimple* ts = dynamic_cast<AliFlowTrackSimple*>(obj);
  if (ts) return PassesCuts(ts);
  return kFALSE; //default when passed a wrong type of object
}

//----------------------------------------------------------------------- 
Bool_t AliFlowTrackSimpleCuts::PassesCuts(const AliFlowTrackSimple *track) const
{
  //simple method to check if the simple track passes the simple cuts
  if(fCutPt) {if (track->Pt() < fPtMin || track->Pt() >= fPtMax ) return kFALSE;}
  if(fCutEta) {if (track->Eta() < fEtaMin || track->Eta() >= fEtaMax ) return kFALSE;}
  if(fCutPhi) {if (track->Phi() < fPhiMin || track->Phi() >= fPhiMax ) return kFALSE;}
  if(fCutCharge) {if (track->Charge() != fCharge) return kFALSE;}
  if(fCutMass) {if (track->Mass() < fMassMin || track->Mass() >= fMassMax ) return kFALSE;}
  //if(fCutPID) {if (track->PID() != fPID) return kFALSE;}
  return kTRUE;
}

//----------------------------------------------------------------------- 
Bool_t AliFlowTrackSimpleCuts::PassesCuts(TParticle* track) const
{
  //simple method to check if the simple track passes the simple cuts
  if(fCutPt)  {if (track->Pt() < fPtMin || track->Pt() >= fPtMax ) return kFALSE;}
  if(fCutEta) {if (track->Eta() < fEtaMin || track->Eta() >= fEtaMax ) return kFALSE;}
  if(fCutPhi) {if (track->Phi() < fPhiMin || track->Phi() >= fPhiMax ) return kFALSE;}
  //if(fCutPID) {if (track->GetPdgCode() != fPID) return kFALSE;}

  //getting the charge from a tparticle is expensive
  //only do it if neccesary
  if (fCutCharge) 
  {
    TParticlePDG* ppdg = track->GetPDG();
    Int_t charge = TMath::Nint(ppdg->Charge()/3.0); //mc particles have charge in units of 1/3e
    return (charge==fCharge);
  }

  if (fCutMass) {
    TParticlePDG* ppdg = track->GetPDG();
    if (ppdg->Mass() < fMassMin || ppdg->Mass() >= fMassMax )
      return kFALSE;
  }

  return kTRUE;
}
