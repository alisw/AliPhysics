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
 
//-------------------------------------------------------------------------
//                      Class AliRsnPID
//                     -------------------
//           Simple collection of reconstructed tracks
//           selected from an ESD event
//           to be used for analysis.
//           .........................................
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#include <TMath.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>

#include "AliLog.h"
#include "AliRsnParticle.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"

#include "AliRsnPID.h"

ClassImp(AliRsnPID)

const char* AliRsnPID::fgkParticleName[AliRsnPID::kSpecies + 1] = 
{
  "electron",
  "muon",
  "pion",
  "kaon",
  "proton",
  "unknown"
};

const Int_t AliRsnPID::fgkParticlePDG[AliRsnPID::kSpecies + 1] = 
{
    11,
    13,
    211,
    321,
    2212,
    0
};

//_____________________________________________________________________________
AliRsnPID::AliRsnPID() :
  fMethod(kRealistic),
  fMaxPt(10.0),
  fMinProb(0.0)
{
//=========================================================
// Default constructor 
// Sets a default setup:
//  - realistic PID
//  - no lower limit for probability
//  - upper limit of 10 GeV for Pt
//  - suitable values for prior probabilities
//=========================================================
    
    fPrior[kElectron] = 0.20;
	fPrior[kMuon    ] = 0.20;
	fPrior[kPion    ] = 0.83;
	fPrior[kKaon    ] = 0.07;
	fPrior[kProton  ] = 0.06;
}

//_____________________________________________________________________________
AliRsnPID::AliRsnPID(const AliRsnPID &event) :
  TObject((TObject)event),
  fMethod(event.fMethod),
  fMaxPt(event.fMaxPt),
  fMinProb(event.fMinProb)
{
//=========================================================
// Copy constructor.
// Implemented to manage the array safely.
//=========================================================
    
    Int_t i;
    for (i = 0; i < kSpecies; i++) fPrior[i] = event.fPrior[i];
}

//_____________________________________________________________________________
AliRsnPID& AliRsnPID::operator=(const AliRsnPID &event)
{
//=========================================================
// Assignment operator.
// Implemented to manage the array safely.
//=========================================================
    
    fMethod = event.fMethod;
    fMaxPt = event.fMaxPt;
    fMinProb = event.fMinProb;
    
    Int_t i;
    for (i = 0; i < kSpecies; i++) fPrior[i] = event.fPrior[i];
    
    // return the newly created object
    return (*this);
}

//_____________________________________________________________________________
AliRsnPID::EType AliRsnPID::InternalType(Int_t pdg)
//=========================================================
// Return the internal enum value corresponding to the PDG
// code passed as argument, if possible.
// Otherwise, returns 'kUnknown' by default.
//=========================================================
{
    EType value;
    Int_t absPDG = TMath::Abs(pdg);
    
	switch (absPDG) {
        case 11:
            value = kElectron;
            break;
        case 13:
            value = kMuon;
            break;
        case 211:
            value = kPion;
            break;
        case 321:
            value = kKaon;
            break;
        case 2212:
            value = kProton;
            break;
        default:
            value = kUnknown;
    }
    return value;
}


//_____________________________________________________________________________
Int_t AliRsnPID::PDGCode(EType type)
{
//=========================================================
// Returns the PDG code of the particle type
// specified as argument (w.r. to the internal enum)
//=========================================================
    
    if (type >= kElectron && type <= kUnknown) {
        return fgkParticlePDG[type];
    }
    else {
        return 0;
    }
}

//_____________________________________________________________________________
const char* AliRsnPID::ParticleName(EType type)
{
//=========================================================
// Returns the name of the particle type
// specified as argument (w.r. to the internal enum)
//=========================================================

    if (type >= kElectron && type <= kUnknown) {
        return fgkParticleName[type];
    }
    else {
        return "unknown";
    }
}

//_____________________________________________________________________________
Double_t AliRsnPID::ParticleMass(EType type)
{
//=========================================================
// Returns the mass corresponding to the particle type
// specified as argument (w.r. to the internal enum)
//=========================================================
    TDatabasePDG *db = TDatabasePDG::Instance();
    Int_t pdg = PDGCode(type);
    return db->GetParticle(pdg)->Mass();
}

//_____________________________________________________________________________
Bool_t AliRsnPID::IdentifyRealistic(AliRsnDaughter *daughter)
{
//=========================================================
// Uses the Bayesian combination of prior probabilities
// with the PID weights of the passed object to compute
// the overall PID probabilities for each particle type.
//
// Once this computation is done, the argument is assigned
// the PID corresponding to the largest probability,
// and its data members are updated accordingly.
// If the track Pt is larger than the cut defined (fMaxPt)
// or the probability is smaller than the cut defined (fMinProb),
// the track is considered unidentified.
//
// If the identification procedure encounters errors, 
// the return value will be "FALSE", otherwise it is "TRUE".
//=========================================================

    // retrieve PID weights from argument
    Int_t i;
    Double_t *prob   = new Double_t[kSpecies];
    Double_t *weight = new Double_t[kSpecies];
    for (i = 0; i < kSpecies; i++) weight[i] = (daughter->PID())[i];
    	
    // step 1 - compute the normalization factor
    Double_t sum = 0.0;
    for (i = 0; i < kSpecies; i++) {
        prob[i] = fPrior[i] * weight[i];
        sum += prob[i];
    }
    if (sum <= 0.0) {
        AliError(Form("Sum of weights = %f < 0.0", sum));
        Unidentify(daughter);
        return kFALSE;
    }
	
	// step 2 - normalize PID weights
    for (i = 0; i < kSpecies; i++) {
        prob[i] /= sum;
    }
	
	// step 3 - finds the maximum probability
    Int_t imax = 0;
    for (i = 1; i < kSpecies; i++) {
        if (prob[i] > prob[imax]) imax = i;
    }
    
    // step 4 - update daughter data members
    //          and check pt & prob
    Double_t pt      = daughter->Pt();
    EType    type    = (EType)imax;
    Double_t maxProb = prob[imax];
    if (pt <= fMaxPt && maxProb >= fMinProb) {
        daughter->SetPIDType(type);
        daughter->SetPIDProb(maxProb);
        daughter->SetM(ParticleMass(type));
    }
    else {
        Unidentify(daughter);
    }
    
    return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnPID::IdentifyPerfect(AliRsnDaughter *daughter)
{
//=========================================================
// Uses the true PDG code to make a perfect identification.
// If the true PDG code does not correspond to any 
// of the expected PID types, gives a warning, sets the 
// PID to 'unknown' and returns kFALSE.
// Otherwise, returns kTRUE.
//=========================================================

    // retrieve true PDG code from argument
    AliRsnParticle *particle = daughter->GetParticle();
    if (!particle) {
        AliWarning("Particle object not initialized: impossible to do perfect PID");
        Unidentify(daughter);
        return kFALSE;
    }
    Int_t pdgCode = particle->PDG();
    EType type = InternalType(pdgCode);
    if (type == kUnknown) {
        //AliWarning(Form("Unrecognized PDG code = %d -- set PID to Unknown", pdgCode));
        Unidentify(daughter);
        return kFALSE;
    }
    else {
        daughter->SetPIDType(type);
        daughter->SetPIDProb(1.0);
        daughter->SetM(ParticleMass(type));
        return kTRUE;
    }
}

//_____________________________________________________________________________
Bool_t AliRsnPID::Unidentify(AliRsnDaughter *daughter)
{
//=========================================================
// Sets the PID to 'unknown' to every track.
//=========================================================

    daughter->SetPIDType(kUnknown);
    return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnPID::Identify(AliRsnDaughter *daughter)
{
//=========================================================
// Recalls one of the above methods, according to the one
// defined in the related data member.
// If the method is not recognized, returns kFALSE and 
// gives an alert. Otherwise, returns kTRUE.
//=========================================================

    switch (fMethod) {
        case kNone:
            Unidentify(daughter);
            return kTRUE;
        case kRealistic:
            IdentifyRealistic(daughter);
            return kTRUE;
        case kPerfect:
            IdentifyPerfect(daughter);
            return kTRUE;
        default:
            AliError(Form("PID method '%d' unrecognized. Nothing done.", fMethod));
            return kFALSE;
    }
}

//_____________________________________________________________________________
Bool_t AliRsnPID::Identify(AliRsnEvent *event)
{
//=========================================================
// Performs identification for all tracks in a given event.
// Returns the logical AND of all PID operations.
//=========================================================

    Bool_t check = kTRUE;
    AliRsnDaughter *daughter = 0;
    TObjArrayIter iter(event->GetAllTracks());
    while ( (daughter = (AliRsnDaughter*)iter.Next()) ) {
        check = check && Identify(daughter);
    }
    event->FillPIDArrays();

    return check;
}


//_____________________________________________________________________________
void AliRsnPID::SetPriorProbability(EType type, Double_t p)
{
//=========================================================
// Sets the prior probability for Realistic PID, for a
// given particle species.
//=========================================================
    
    if (type >= kElectron && type < kSpecies) {
        fPrior[type] = p;
    }
}
