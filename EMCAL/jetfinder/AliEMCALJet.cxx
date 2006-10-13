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

//*-- Author: Andreas Morsch (CERN)

#include "AliEMCALJet.h"
#include "Ecommon.h"

ClassImp(AliEMCALJet)

//____________________________________________________________________________
AliEMCALJet::AliEMCALJet()
  : fEnergy(0.), fEMCALEnergy(0.),
    fEMCALEnergyBGSub(0), fTrackEnergy(0.),
    fTrackEnergyPtCut(0.), fHCEnergy(0.),
    fIsWeightedEnergy(kFALSE), fEta(0.),
    fPhi(0.), fNt(0), fPtT(0), fEtaT(0),
    fPhiT(0), fPdgT(0)
{
// Default constructor
}

AliEMCALJet::AliEMCALJet(Float_t energy, Float_t phi, Float_t eta)
  : fEnergy(energy), fEMCALEnergy(0.),
    fEMCALEnergyBGSub(0), fTrackEnergy(0.),
    fTrackEnergyPtCut(0.), fHCEnergy(0.),
    fIsWeightedEnergy(kFALSE), fEta(eta),
    fPhi(phi), fNt(0), fPtT(0), fEtaT(0),
    fPhiT(0), fPdgT(0)
{
// Constructor
}

AliEMCALJet::AliEMCALJet(const AliEMCALJet& jet) 
  : TObject(jet), fEnergy(jet.fEnergy), fEMCALEnergy(jet.fEMCALEnergy),
    fEMCALEnergyBGSub(jet.fEMCALEnergyBGSub), 
    fTrackEnergy(jet.fTrackEnergy),
    fTrackEnergyPtCut(jet.fTrackEnergyPtCut), 
    fHCEnergy(jet.fHCEnergy),
    fIsWeightedEnergy(jet.fIsWeightedEnergy), 
    fEta(jet.fEta),fPhi(jet.fPhi), fNt(jet.fNt), 
    fPtT(jet.fPtT), fEtaT(jet.fEtaT),
    fPhiT(jet.fPhiT), fPdgT(jet.fPdgT)
{
// Copy Constructor
}

//____________________________________________________________________________
AliEMCALJet::~AliEMCALJet()
{
// Destructor
}


void AliEMCALJet::SetTrackList(Int_t n, Float_t* pt, Float_t* eta, Float_t* phi, Int_t* pdg)
{
//
// 
    fNt = n;
    for (Int_t i = 0; i < n; i++) {
	fPtT [i]  = pt [i];
	fEtaT[i]  = eta[i];
	fPhiT[i]  = phi[i];
	fPdgT[i]  = pdg[i];
    }
}



Int_t AliEMCALJet::TrackList(Float_t* pt, Float_t* eta, Float_t* phi, Int_t* pdg) const
{
//
// 
    for (Int_t i = 0; i < fNt; i++) {
	pt [i] = fPtT [i];
	eta[i] = fEtaT[i];
	phi[i] = fPhiT[i];
	pdg[i] = fPdgT[i];
    }
    return fNt;
}







