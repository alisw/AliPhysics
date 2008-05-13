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

// ***********************
// *** AliRsnDaughter ****
// ***********************
// 
// Light-weight 'track' object into an internal format used
// for further steps of resonance analysis.
// Provides converters from all kinds of input track type
// (ESD, AOD and MC).
// Contains also a facility to compute invariant mass of a pair.
//
// author: A. Pulvirenti --- email: alberto.pulvirenti@ct.infn.it

#include <Riostream.h>

#include <TParticle.h>
#include <TString.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliMCParticle.h"

#include "AliRsnPID.h"
#include "AliRsnParticle.h"
#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter() :
  AliVParticle(),
  fIndex(-1),
  fLabel(-1),
  fCharge(0),
  fFlags(0),
  fPIDType(AliRsnPID::kUnknown),
  fPIDProb(0.0),
  fMass(0.0),
  fParticle(0x0)
{
//=========================================================
// Default constructor.
// Initializes all data-members with meaningless values.
//=========================================================

    Int_t i;
    for (i = 0; i < AliRsnPID::kSpecies; i++) {
        if (i < 3) {
            fP[i] = 0.0;
            fV[i] = 0.0;
        }
        fPIDWeight[i] = 0.0;
    }
}

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(const AliRsnDaughter &copy) :
  AliVParticle(copy),
  fIndex(copy.fIndex),
  fLabel(copy.fLabel),
  fCharge(copy.fCharge),
  fFlags(copy.fFlags),
  fPIDType(copy.fPIDType),
  fPIDProb(copy.fPIDProb),
  fMass(copy.fMass),
  fParticle(0x0)
{
//=========================================================
// Copy constructor.
//=========================================================

    Int_t i;
    for (i = 0; i < AliRsnPID::kSpecies; i++) {
        if (i < 3) {
            fP[i] = copy.fP[i];
            fV[i] = copy.fV[i];
        }
        fPIDWeight[i] = copy.fPIDWeight[i];
    }

    // initialize particle object 
    // only if it is present in the template object
    if (copy.fParticle) fParticle = new AliRsnParticle(*(copy.fParticle));
}

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(AliESDtrack *track) :
  AliVParticle(),
  fIndex(-1),
  fLabel(-1),
  fCharge(0),
  fFlags(0),
  fPIDType(AliRsnPID::kUnknown),
  fPIDProb(0.0),
  fMass(0.0),
  fParticle(0x0)
{
//=========================================================
// Constructor to get data from an ESD track.
//=========================================================

    Adopt(track);
}

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(AliAODTrack *track) :
  AliVParticle(),
  fIndex(-1),
  fLabel(-1),
  fCharge(0),
  fFlags(0),
  fPIDType(AliRsnPID::kUnknown),
  fPIDProb(0.0),
  fMass(0.0),
  fParticle(0x0)
{
//=========================================================
// Constructor to get data from an AOD track.
//=========================================================

    Adopt(track);
}

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(AliMCParticle *track) :
  AliVParticle(),
  fIndex(-1),
  fLabel(-1),
  fCharge(0),
  fFlags(0),
  fPIDType(AliRsnPID::kUnknown),
  fPIDProb(0.0),
  fMass(0.0),
  fParticle(0x0)
{
//=========================================================
// Constructor to get data from an MC track.
//=========================================================

    Adopt(track);
}

//_____________________________________________________________________________
AliRsnDaughter& AliRsnDaughter::operator=(const AliRsnDaughter &copy)
{
//=========================================================
// Assignment operator.
// Works like the copy constructor and returns a reference
// to the initialized object for which it is called.
//=========================================================
    
    fIndex  = copy.fIndex;
    fLabel  = copy.fLabel;
    fCharge = copy.fCharge;
    fFlags  = copy.fFlags;

    Int_t i;
    for (i = 0; i < AliRsnPID::kSpecies; i++) {
        if (i < 3) {
            fP[i] = copy.fP[i];
            fV[i] = copy.fV[i];
        }
        fPIDWeight[i] = copy.fPIDWeight[i];
    }
    
    fPIDType = copy.fPIDType;
    fPIDProb = copy.fPIDProb;
    fMass    = copy.fMass;
    
    // initialize particle object 
    // only if it is present in the template object;
    // otherwise, it is just cleared and not replaced with anything
    if (fParticle) {
        delete fParticle;
        fParticle = 0x0;
    }
    if (copy.fParticle) fParticle = new AliRsnParticle(*(copy.fParticle));

    return (*this);
}

//_____________________________________________________________________________
AliRsnDaughter::~AliRsnDaughter()
{
//=========================================================
// Destructor
//=========================================================

    if (fParticle) {
        delete fParticle;
        fParticle = 0;
    }
}

//_____________________________________________________________________________
void AliRsnDaughter::SetPIDWeight(Int_t i, Double_t value)
{
//=========================================================
// I the argument 'i' is in the correct range,
// sets the i-th PID weight to 'value'
//=========================================================

    if (i >= 0 && i < AliRsnPID::kSpecies) fPIDWeight[i] = value;
    else {
        AliError(Form("Cannot set a weight related to slot %d", i));
    }
}


//_____________________________________________________________________________
void AliRsnDaughter::SetPIDWeights(const Double_t *pid)
{
//=========================================================
// Sets ALL PID weights at once.
// The passed array is supposed to have at least as many
// slots as the number of allowed particle species.
//=========================================================

   Int_t i;
   for (i = 0; i < AliRsnPID::kSpecies; i++) fPIDWeight[i] = pid[i];
}


//_____________________________________________________________________________
Bool_t AliRsnDaughter::Adopt(AliESDtrack* esdTrack)
{
//=========================================================
// Copies data from an AliESDtrack into "this":
//  - charge sign
//  - momentum
//  - point of closest approach to primary vertex
//  - ESD pid weights
// In case of errors returns kFALSE, otherwise kTRUE.
//=========================================================

    if (!esdTrack) {
        AliError("Passed NULL object: nothing done");
        return kFALSE;
    }
	
    // copy momentum and vertex
    esdTrack->GetPxPyPz(fP);
    esdTrack->GetXYZ(fV);

    // copy PID weights
    Int_t    i;
    Double_t pid[5];
    esdTrack->GetESDpid(pid);
    for (i = 0; i < 5; i++) fPIDWeight[i] = pid[i];
    
    // copy flags
    fFlags = esdTrack->GetStatus();
    
    // copy charge sign
    fCharge = (Short_t)esdTrack->Charge();
 
    return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliRsnDaughter::Adopt(AliAODTrack* aodTrack)
{
//=========================================================
// Copies data from an AliAODtrack into "this":
//  - charge sign
//  - momentum
//  - point of closest approach to primary vertex
//  - ESD pid weights
// In case of errors returns kFALSE, otherwise kTRUE.
//=========================================================

    if (!aodTrack) {
        AliError("Passed NULL object: nothing done");
        return kFALSE;
    }
	
    // copy momentum  and vertex
    fP[0] = aodTrack->Px();
    fP[1] = aodTrack->Py();
    fP[2] = aodTrack->Pz();
    fV[0] = aodTrack->Xv();
    fV[1] = aodTrack->Yv();
    fV[2] = aodTrack->Zv();
    
    // copy PID weights
    Int_t i;
    for (i = 0; i < 5; i++) fPIDWeight[i] = aodTrack->PID()[i];
    
    // copy sign
    fCharge = aodTrack->Charge();

    return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliRsnDaughter::Adopt(AliMCParticle *mcParticle)
{
//=========================================================
// Copies data from a MCParticle into "this":
//  - charge sign
//  - momentum
//  - point of closest approach to primary vertex
//  - ESD pid weights
//  - true PDG code
//  - mother
// In case of errors returns kFALSE, otherwise kTRUE.
//=========================================================

    if (!mcParticle) {
        AliError("Passed NULL object: nothing done");
        return kFALSE;
    }
	
	// retrieve the TParticle object from the argument
	TParticle *particle = mcParticle->Particle();
	if (!particle) {
	   AliError("AliMCParticle::Particle() returned NULL");
	   return kFALSE;
    }
	
    // copy momentum  and vertex
    fP[0] = particle->Px();
    fP[1] = particle->Py();
    fP[2] = particle->Pz();
    fV[0] = particle->Vx();
    fV[1] = particle->Vy();
    fV[2] = particle->Vz();
    
    // recognize charge sign from PDG code sign
    Int_t pdg = particle->GetPdgCode();
    Int_t absPDG = TMath::Abs(pdg);
    if (absPDG <= 15) {
        if (pdg > 0) fCharge = -1; else fCharge = 1;
    }
    else if (absPDG < 3000) {
        if (pdg > 0) fCharge = 1; else fCharge = -1;
    }
    else {
        fCharge = 0;
        return kFALSE;
    }
    
    // identify track perfectly using PDG code
    fPIDType = AliRsnPID::InternalType(pdg);
    fPIDProb = 1.0;
    fMass = AliRsnPID::ParticleMass(fPIDType);
    
    // flags and PID weights make no sense with MC tracks
    fFlags = 0;
    for (pdg = 0; pdg < AliRsnPID::kSpecies; pdg++) fPIDWeight[pdg] = 0.0;
    fPIDWeight[fPIDType] = 1.0;
    
    // copy other MC info (mother PDG code cannot be retrieved here)
    InitParticle(particle);

    return kTRUE;
}

//_____________________________________________________________________________
void AliRsnDaughter::Print(Option_t *option) const
{
//=========================================================
// Prints the values of data members, using the options:
// - P --> momentum
// - V --> DCA vertex
// - C --> electric charge
// - F --> flags
// - I --> identification (PID, probability and mass)
// - W --> PID weights
// - M --> Montecarlo (from AliRsnParticle)
// - ALL --> All oprions switched on
//
// Index and label are printed by default.
//=========================================================

    TString opt(option);
    opt.ToUpper();
    
    cout << ".......Index            : " << fIndex << endl;
    cout << ".......Label            : " << fLabel << endl;
    
    if (opt.Contains("P") || opt.Contains("ALL")) {
        cout << ".......Px, Py, Pz, Pt   : " << Px() << ' ' << Py() << ' ' << Pz() << ' ' << Pt() << endl;
    }
    if (opt.Contains("V") || opt.Contains("ALL")) {
        cout << ".......Vx, Vy, Vz       : " << Xv() << ' ' << Yv() << ' ' << Zv() << endl;
    }
    if (opt.Contains("C") || opt.Contains("ALL")) {
        cout << ".......Charge           : " << fCharge << endl;
    }
    if (opt.Contains("F") || opt.Contains("ALL")) {
        cout << ".......Flags            : " << fFlags << endl;
    }
    if (opt.Contains("I") || opt.Contains("ALL")) {
        cout << ".......PID              : " << AliRsnPID::ParticleName(fPIDType) << endl;
        cout << ".......PID probability  : " << fPIDProb << endl;
        cout << ".......Mass             : " << fMass << endl;
    }
    if (opt.Contains("W") || opt.Contains("ALL")) {
        cout << ".......Weights          : ";
        Int_t i;
        for (i = 0; i < AliRsnPID::kSpecies; i++) cout << fPIDWeight[i] << ' ';
        cout << endl;
    }
    if (opt.Contains("M") || opt.Contains("ALL")) {
        if (fParticle) {
            cout << ".......PDG code         : " << fParticle->PDG() << endl;
            cout << ".......Mother (label)   : " << fParticle->Mother() << endl;
            cout << ".......Mother (PDG code): " << fParticle->MotherPDG() << endl;
        }
        else {
            cout << ".......MC info not present" << endl;
        }
    }
}


//_____________________________________________________________________________
AliRsnDaughter AliRsnDaughter::Sum(AliRsnDaughter t1, AliRsnDaughter t2)
{
//=========================================================
// Creates an AliRsnDaughter object composing momenta
// of the particles passed as arguments.
// The vector momentum of the returned object is 
// the sum of the ones of the arguments.
// The mass of the object is set to the invariant mass
// computed from the relativistic sum of the 4-momenta
// of the two arguments.
// Finally, if both arguments have the same 'fMother' 
// (when this MC information is available), it is set
// also as the 'fMother' of the returned object, together
// with its PDG code.
//=========================================================

    // create a 'default' new AliRsnDaughter
    AliRsnDaughter out;
    
    // define the momentum of the sum as the 
    // relativistic sum of 4-momenta of arguments
    Double_t etot  = t1.E()  + t2.E();
    Double_t pxTot = t1.Px() + t2.Px();
    Double_t pyTot = t1.Py() + t2.Py();
    Double_t pzTot = t1.Pz() + t2.Pz();
    Double_t mass  = TMath::Sqrt(etot*etot - pxTot*pxTot - pyTot*pyTot - pzTot*pzTot);
    out.SetM(mass);
    out.SetP(pxTot, pyTot, pzTot);
    
    // if both arguments have MC info, it is used here
    // to check if they (in the truth MC) are daughters
    // of the same resonance and, if this is the case,
    // details of that resonance are set into the 
    // 'fParticle' object of the sum
    Int_t mum1, mum2;
    if (t1.GetParticle() && t2.GetParticle()) {
        out.InitParticle();
        mum1 = t1.GetParticle()->Mother();
        mum2 = t2.GetParticle()->Mother();
        if (mum1 == mum2 && mum1 >= 0) {
            out.SetLabel(mum1);
            out.GetParticle()->SetPDG(t1.GetParticle()->MotherPDG());
        }
    }

    return out;
}


//_____________________________________________________________________________
void AliRsnDaughter::InitParticle()
{
//=========================================================
// Initializes the particle object with default constructor.
//=========================================================

    fParticle = new AliRsnParticle;
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::InitParticle(TParticle *particle)
{
//=========================================================
// Copies data from an MC particle into the object which
// contains all MC details taken from kinematics info.
// If requested by second argument, momentum and vertex
// of the Particle are copied into the 'fP' and 'fV'
// data members, to simulate a perfect reconstruction.
// If something goes wrong, returns kFALSE, 
// otherwise returns kTRUE.
//=========================================================

    // retrieve the TParticle object pointed by this MC track
    if (!particle) {
        AliError("Passed NULL particle object");
        return kFALSE;
    }
	
    // initialize object if not initialized yet
    if (fParticle) delete fParticle;
    fParticle = new AliRsnParticle;
    fParticle->Adopt(particle);

    return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliRsnDaughter::InitParticle(AliMCParticle *mcParticle)
{
//=========================================================
// Copies data from an MC particle into the object which
// contains all MC details taken from kinematics info.
// If requested by second argument, momentum and vertex
// of the Particle are copied into the 'fP' and 'fV'
// data members, to simulate a perfect reconstruction.
// If something goes wrong, returns kFALSE, 
// otherwise returns kTRUE.
//=========================================================

    // retrieve the TParticle object pointed by this MC track
    TParticle *particle = mcParticle->Particle();
    return InitParticle(particle);
}
