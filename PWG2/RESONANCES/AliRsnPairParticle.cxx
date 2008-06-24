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

//
// Class AliRsnPairParticle
//
// Implementation of a pair of tracks, for several purposes
// - computing the total 4-momentum & inv. mass for output histos filling
// - evaluating cut checks on the pair of particles
//
// author: Martin Vala (martin.vala@cern.ch)
// revised by: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliLog.h"

#include "AliRsnPairParticle.h"

ClassImp (AliRsnPairParticle)

//_____________________________________________________________________________
AliRsnPairParticle::AliRsnPairParticle() :
  fMass(0.0)
{
//
// Constructor.
// Initializes all variables to meaningless values.
//
    Int_t i, j;
    for (i = 0; i < 3; i++) {
        fPTot[i] = 0.0;
        fPTotMC[i] = 0.0;
        if (i < 2) {
            fMotherLabel[i] = -1;
            fMotherPDG[i] = 0;
            fDaughter[i] = 0x0;
        }
        for (j = 0; j < 2; j++) {
            fPTrack[j][i] = 0.0;
            fPTrackMC[j][i] = 0.0;
        }
    }
}

//_____________________________________________________________________________
AliRsnPairParticle::AliRsnPairParticle(const AliRsnPairParticle &obj) :
  TObject(obj),
  fMass(obj.fMass)
{
//
// Copy constructor.
// Initializes all variables to copy values.
// Does not duplicate pointers.
//
    Int_t i, j;
    for (i = 0; i < 3; i++) {
        fPTot[i] = obj.fPTot[i];
        fPTotMC[i] = obj.fPTotMC[i];
        if (i < 2) {
            fMotherLabel[i] = obj.fMotherLabel[i];
            fMotherPDG[i] = obj.fMotherPDG[i];
            fDaughter[i] = obj.fDaughter[i];
        }
        for (j = 0; j < 2; j++) {
            fPTrack[j][i] = obj.fPTrack[j][i];
            fPTrackMC[j][i] = obj.fPTrackMC[j][i];
        }
    }
}

//_____________________________________________________________________________
AliRsnPairParticle& AliRsnPairParticle::operator=(const AliRsnPairParticle &obj)
{
//
// Assignment operator.
// Initializes all variables to copy values.
// Does not duplicate pointers.
//

    fMass = obj.fMass;

    Int_t i, j;
    for (i = 0; i < 3; i++) {
        fPTot[i] = obj.fPTot[i];
        fPTotMC[i] = obj.fPTotMC[i];
        if (i < 2) {
            fMotherLabel[i] = obj.fMotherLabel[i];
            fMotherPDG[i] = obj.fMotherPDG[i];
            fDaughter[i] = obj.fDaughter[i];
        }
        for (j = 0; j < 2; j++) {
            fPTrack[j][i] = obj.fPTrack[j][i];
            fPTrackMC[j][i] = obj.fPTrackMC[j][i];
        }
    }

    return (*this);
}

//_____________________________________________________________________________
AliRsnPairParticle::~AliRsnPairParticle()
{
//
// Desctructor.
// Does nothing.
//
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetInvMass(Double_t m1, Double_t m2)
{
//
// Compute invariant mass using REC values.
// The masses in argument have default value = -1.0;
// when this is used, the mass for computation is taken from daughter object.
// Otherwise, the passed values are used.
//
    Double_t etot = 0.0;

    // energy of particle 1
    if (m1 > 0.0) etot = GetDaughterEnergy(0, m1);
    else etot = fDaughter[0]->E();

    // energy of particle 2
    if (m2 > 0.0) etot += GetDaughterEnergy(1, m2);
    else etot += fDaughter[1]->E();

    // total momentum square module
    Double_t p2Tot = fPTot[0]*fPTot[0] + fPTot[1]*fPTot[1] + fPTot[2]*fPTot[2];

    // invariant mass returned
    return  TMath::Sqrt (etot * etot - p2Tot);
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetInvMassMC(Double_t m1, Double_t m2)
{
//
// Compute invariant mass using MC values.
// The masses in argument have default value = -1.0;
// when a negative mass value is passed as one of the arguments,
// the energy of tracks for computation is taken from daughter object.
// Otherwise, the passed mass values are used to re-compute the track energy.
//
    Double_t etot = 0.0;

    // energy of particle 1
    if (m1 > 0.0) etot = GetDaughterEnergyMC(0, m1);
    else etot = fDaughter[0]->GetMCInfo()->E();

    // energy of particle 2
    if (m2 > 0.0) etot += GetDaughterEnergyMC(1, m2);
    else etot += fDaughter[1]->GetMCInfo()->E();

    // total momentum square module
    Double_t p2Tot = fPTotMC[0]*fPTotMC[0] + fPTotMC[1]*fPTotMC[1] +fPTotMC[2]*fPTotMC[2];

    // invariant mass returned
    return  TMath::Sqrt (etot * etot - p2Tot);
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetDaughterEnergy(const Int_t &index, const Double_t &mass)
{
//
// Compute track energy from REC momentum using a passed mass value.
// The index argument refers to the used track among the two of the pair.
//

    if (mass > 0 && index >= 0 && index < 2) {
        AliRsnDaughter temp(*fDaughter[index]);
        temp.SetM(mass);
        return temp.E();
    }

    AliWarning("Negative mass or wrong index passed. Returning 0");
    return 0.0;
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetDaughterEnergyMC(const Int_t &index, const Double_t &mass)
{
//
// Compute track energy from MC momentum using a passed mass value.
// The index argument refers to the used track among the two of the pair.
//
    Int_t i;
    Double_t p2Tot = 0.0;
    if (mass > 0 && index >= 0 && index < 2) {
        for (i = 0; i < 3; i++) p2Tot += fPTrackMC[index][i] * fPTrackMC[index][i];
        return TMath::Sqrt(mass*mass + p2Tot);
    }

    AliWarning("Negative mass or wrong index passed. Returning 0");
    return 0.0;
}

//_____________________________________________________________________________
Int_t AliRsnPairParticle::GetPDG(const Int_t &index)
{
//
// Return PDG code of one track in the pair
//

    if (index < 0 || index > 1) {
        AliError("Index out of range");
        return 0;
    }
    if (!fDaughter[index]->GetMCInfo()) {
        AliError(Form("MCInfo not initialized for track %d", index));
        return 0;
    }
    return fDaughter[index]->GetMCInfo()->PDG();
}

//_____________________________________________________________________________
Bool_t AliRsnPairParticle::IsTruePair(Int_t refPDG)
{
//
// Checks if the tweo tracks in the pair come from the same resonance.
// This can be known if MC info is present, looking at the GEANT label of mother
// (which should be the same).
// If the argument is 0, the answer is kTRUE whenever the labels of mothers of the
// two tracks is the same. When the argument is not zero, the answer is kTRUE only
// if the mother is the same and its PDG code is equal to the argument.

    // if MC info is not available, the pairs is not true by default
    if (!fDaughter[0]->GetMCInfo() || !fDaughter[1]->GetMCInfo()) {
        return kFALSE;
    }

    // check that labels are the same
    if (fDaughter[0]->GetMCInfo()->Mother() != fDaughter[1]->GetMCInfo()->Mother()) {
        return kFALSE;
    }

    // if we reach this point, the two tracks have the same mother
    // let's check now the PDG code of this common mother
    Int_t motherPDG = TMath::Abs(fDaughter[0]->GetMCInfo()->MotherPDG());
    if (refPDG == 0) return kTRUE;
    else return (motherPDG == refPDG);
}

//_____________________________________________________________________________
void AliRsnPairParticle::SetPair(AliRsnDaughter *daughter1, AliRsnDaughter *daughter2)
{
//
// Main object filling method.
// Accepts two AliRsnDaughter's which are the two tracks in the pair,
// fills all data-members which contain their momenta & info,
// and computes the total momentum for REC data and MC if available
//

    Int_t i;

    fDaughter[0] = daughter1;
    fDaughter[1] = daughter2;

    // copy MC info (if available)
    if (fDaughter[0]->GetMCInfo() && fDaughter[1]->GetMCInfo()) {
        for (i = 0; i < 2; i++) {
            fPTrackMC[i][0] = fDaughter[i]->GetMCInfo()->Px();
            fPTrackMC[i][1] = fDaughter[i]->GetMCInfo()->Py();
            fPTrackMC[i][2] = fDaughter[i]->GetMCInfo()->Pz();
            fMotherPDG[i] = fDaughter[i]->GetMCInfo()->MotherPDG();
        }
        for (i = 0; i < 3; i++) fPTotMC[i] = fPTrackMC[0][i] + fPTrackMC[1][i];
    }

    // copy reconstructed info (always available)
    for (i = 0; i < 2; i++) {
        fPTrack[i][0] = fDaughter[i]->Px();
        fPTrack[i][1] = fDaughter[i]->Py();
        fPTrack[i][2] = fDaughter[i]->Pz();
    }
    for (i = 0; i < 3; i++) fPTot[i] = fPTrack[0][i] + fPTrack[1][i];
}

//_____________________________________________________________________________
void AliRsnPairParticle::PrintInfo (const Option_t *option)
{
//
// Print some info of the pair
//

    TString s(option);

    AliInfo ( "======== BEGIN PAIR INFO ===========" );
    if (s.Contains("p")) {
        AliInfo (Form("Px1 = %.6f --- Py1 = %.6f --- Pz1 = %.6f", fPTrack[0][0], fPTrack[0][1], fPTrack[0][2]));
        AliInfo (Form("Px2 = %.6f --- Py2 = %.6f --- Pz2 = %.6f", fPTrack[1][0], fPTrack[1][1], fPTrack[1][2]));
    }
    if (s.Contains("t")) {
        AliInfo (Form("PDG1 = %d --- PDG2 = %d", GetPDG(0), GetPDG(1)));
        AliInfo (Form("type1 = %d --- type2 = %d", GetType(0), GetType(1)));
        AliInfo (Form("label1 = %d --- label2 = %d", GetLabel(0), GetLabel(1)));
    }
    AliInfo ( "====== END PAIR INFO =========" );
}
