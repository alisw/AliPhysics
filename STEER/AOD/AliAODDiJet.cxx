/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     AOD class for di-jets
//     The present version is for test purposes only
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TLorentzVector.h>
#include <TProcessID.h>
#include "AliAODDiJet.h"
#include "AliAODJet.h"

ClassImp(AliAODDiJet)


//______________________________________________________________________________
AliAODDiJet::AliAODDiJet() :
    AliAODJet(),
    fJetR(0),
    fJet1(0),
    fJet2(0)
{
  // constructor
}

AliAODDiJet::AliAODDiJet(Double_t px, Double_t py, Double_t pz, Double_t e):
    AliAODJet(px, py, pz, e), 
    fJetR(0),
    fJet1(0),
    fJet2(0)
{
  // another constructor
}

AliAODDiJet::AliAODDiJet(TLorentzVector & p):
    AliAODJet(p),
    fJetR(0),
    fJet1(0),
    fJet2(0)
{
  // constructor
}


//______________________________________________________________________________
AliAODDiJet::~AliAODDiJet() 
{
  // destructor
    if (fJetR) fJetR->Delete();
}

void AliAODDiJet::SetJetRefs(AliAODJet* jet1, AliAODJet* jet2) 
{
// Set references to the two jets
    if (fJetR) fJetR->Delete();
    fJetR = new TRefArray(TProcessID::GetProcessWithUID( jet1 ));
    fJetR->Add(jet1);
    fJetR->Add(jet2);
    fJet1 = jet1;
    fJet2 = jet2;
}


//______________________________________________________________________________
Float_t AliAODDiJet::DeltaPhi()
{
  // phi distance between the two jets
  // the result is in the interval [0,pi]

    Float_t phi1 = Jet(0)->Phi();
    Float_t phi2 = Jet(1)->Phi();
    Float_t dphi = TMath::Abs(phi1 - phi2);
    dphi = (dphi > TMath::Pi()) ? (2*TMath::Pi()-dphi) : dphi;
    return dphi;
}


//______________________________________________________________________________
Float_t AliAODDiJet::PhiImbalance()
{
  // phi imbalance wrt back-to-back
  // the result is in the interval [-pi,pi]

    Float_t phi1 = Jet(0)->Phi();
    Float_t phi2 = Jet(1)->Phi();
    Float_t dphi = phi1 - phi2;
    dphi = dphi - TMath::Sign((Float_t)TMath::Pi(),dphi);
    return dphi;
}
