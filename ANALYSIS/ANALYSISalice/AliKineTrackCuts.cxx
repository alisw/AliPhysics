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

#include <TObject.h>
#include <TMath.h>
#include <TParticle.h>

#include "AliKineTrackCuts.h"

//
//  Class for simple Kinematic cuts on
//  particles (tracks) from Kinematic stack (TParticle)
//  MC Simulation
//


//____________________________________________________________________
ClassImp(AliKineTrackCuts)

//____________________________________________________________________
AliKineTrackCuts::AliKineTrackCuts(const Char_t* name, const Char_t* title) : 
  AliAnalysisCuts(name,title),
  fOnlyFinalParticles(kFALSE),
  fOnlyPrimary(kFALSE), 
  fPMin(0),
  fPMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fRapMin(0),
  fRapMax(0)
{
  //
  // constructor
  //
  // setting default cuts
  SetPRange();
  SetPtRange();
  SetPxRange();
  SetPyRange();
  SetPzRange();
  SetEtaRange();
  SetRapRange();
}



//____________________________________________________________________
Bool_t  AliKineTrackCuts::IsSelected(TObject* obj)
{

  TParticle * part = (TParticle *)obj;
  
  // only final particles
  if( fOnlyFinalParticles && part->GetStatusCode() !=1 ) return kFALSE;
  if( fOnlyPrimary && part->IsPrimary() !=1 ) return kFALSE;
  
  // getting the kinematic variables of the track
  Float_t momentum = part->P();
  Float_t pt       = part->Pt();
  Float_t energy   = part->Energy();

  //y-eta related calculations
  Float_t eta = part->Eta();
  Float_t y   = -100.;
  if((energy != TMath::Abs(part->Pz()))&&(momentum != 0))
    y = 0.5*TMath::Log((energy + part->Pz())/(energy - part->Pz()));

  if((momentum < fPMin) || (momentum > fPMax)) return kFALSE;
  if((pt < fPtMin) || (pt > fPtMax)) return kFALSE;
  if((part->Px() < fPxMin) || (part->Px() > fPxMax)) return kFALSE;
  if((part->Py() < fPyMin) || (part->Py() > fPyMax)) return kFALSE;
  if((part->Pz() < fPzMin) || (part->Pz() > fPzMax)) return kFALSE;
  if((eta < fEtaMin) || (eta > fEtaMax)) return kFALSE;
  if((y < fRapMin) || (y > fRapMax)) return kFALSE;

  return kTRUE;
}









