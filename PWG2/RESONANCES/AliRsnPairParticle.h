/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
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

#ifndef ALIRSNPAIRPARTICLE_H
#define ALIRSNPAIRPARTICLE_H

#include <TMath.h>

#include "AliRsnDaughter.h"
#include "AliRsnMCInfo.h"

class AliRsnPairParticle : public TObject
{
public:

    AliRsnPairParticle();
    AliRsnPairParticle(const AliRsnPairParticle &obj);
    AliRsnPairParticle& operator=(const AliRsnPairParticle &obj);
    virtual ~AliRsnPairParticle();

    Double_t   GetInvMass (Double_t m1 = -1.0, Double_t m2 = -1.0);
    Double_t   GetInvMassMC (Double_t m1 = -1.0, Double_t m2 = -1.0);
    Double_t   GetDaughterEnergy(const Int_t &index, const Double_t &mass);
    Double_t   GetDaughterEnergyMC(const Int_t &index, const Double_t &mass);

    Double_t   GetP2() const {return (fPTot[0]*fPTot[0] + fPTot[1]*fPTot[1] + fPTot[2]*fPTot[2]);}
    Double_t   GetPt2() const {return (fPTot[0]*fPTot[0] + fPTot[1]*fPTot[1]);}
    Double_t   GetP() const {return TMath::Sqrt(GetP2());}
    Double_t   GetP(Int_t index) const {return fPTot[index];}
    Double_t   GetPt() const {return TMath::Sqrt(GetPt2());}

    Double_t   GetP2MC() const {return (fPTotMC[0]*fPTotMC[0] + fPTotMC[1]*fPTotMC[1] + fPTotMC[2]*fPTotMC[2]);}
    Double_t   GetPt2MC() const {return (fPTotMC[0]*fPTotMC[0] + fPTotMC[1]*fPTotMC[1]);}
    Double_t   GetPMC() const {return TMath::Sqrt(GetP2MC());}
    Double_t   GetPMC(Int_t index) const {return fPTotMC[index];}
    Double_t   GetPtMC() const {return TMath::Sqrt(GetPt2MC());}

    Int_t      GetPDG(const Int_t &index);
    Int_t      GetType(const Int_t &index) const {return (Int_t) fDaughter[index]->PIDType();}
    Int_t      GetLabel(const Int_t &index) const {return (Int_t) fDaughter[index]->Label();}
    Int_t      GetIndex(const Int_t &index) const {return (Int_t) fDaughter[index]->Index();}
    Double_t   GetMass() const {return fMass;}
    AliRsnDaughter* GetDaughter(const Int_t &index) const {return fDaughter[index];}

    Bool_t     IsPDGEqual() {return abs(GetPDG(0)) == abs(GetPDG(1));}
    Bool_t     IsTypeEqual() {return GetType(0) == GetType(1);}
    Bool_t     IsLabelEqual() {return abs(GetLabel(0)) == abs(GetLabel(1));}
    Bool_t     IsIndexEqual() {return GetIndex(0) == GetIndex(1);}
    Bool_t     IsTruePair(Int_t refPDG = 0);

    void       SetMass(Double_t mass) {fMass = mass;}
    void       SetPair(AliRsnDaughter *daughter1, AliRsnDaughter *daughter2);
    void       PrintInfo (const Option_t *option = "");

private:

    Double_t         fPTot[3];          // total momentum computed with rec. values
    Double_t         fPTotMC[3];        // total momentum computed with MC values
    Double_t         fPTrack[2][3];     // rec. momentum of single tracks
    Double_t         fPTrackMC[2][3];   // MC momentum of single tracks

    Double_t         fMass;             // mass hypothesis for resonance

    Int_t            fMotherLabel[2];   // GEANT label of tracks
    Int_t            fMotherPDG[2];     // PDG code of mother of tracks

    AliRsnDaughter  *fDaughter[2];      // elements of the pair

    ClassDef (AliRsnPairParticle,1)
};

#endif
