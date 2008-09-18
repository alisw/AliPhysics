//
// Class AliRsnPairDef
//
// Defines a decay channel for a resonance,
// resulting in a specified PDG code for the mother,
// and the particle type for the daughters, defined
// according to the internal PID format of the package
//
// author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNPAIRDEF_H
#define ALIRSNPAIRDEF_H

#include <TString.h>

#include "AliRsnPID.h"

class AliRsnDaughter;

class AliRsnPairDef : public TObject
{
  public:

    AliRsnPairDef();
    AliRsnPairDef(Char_t sign1, AliRsnPID::EType type1,
                  Char_t sign2, AliRsnPID::EType type2, Int_t motherPDG = 0);
    AliRsnPairDef(AliRsnPID::EType type1, Char_t sign1,
                  AliRsnPID::EType type2, Char_t sign2, Int_t motherPDG = 0);
    AliRsnPairDef(const AliRsnPairDef &copy);
    const AliRsnPairDef& operator= (const AliRsnPairDef &copy);
    virtual ~AliRsnPairDef() { }

    // getters
    Char_t           GetCharge(Int_t i) const {if (i>=0&&i<2) return fCharge[i]; else return 0;}
    AliRsnPID::EType GetType(Int_t i) const {if (i>=0&&i<2) return fType[i]; else return AliRsnPID::kUnknown;}
    Double_t         GetMass(Int_t i) const {if (i>=0&&i<2) return fMass[i]; else return 0.0;}
    Int_t            GetMotherPDG() const {return fMotherPDG;}
    TString          GetPairName();

    // setters
    Bool_t SetPairElement(Int_t i, Char_t charge, AliRsnPID::EType pid);
    Bool_t SetPair(Char_t ch1, AliRsnPID::EType pid1, Char_t ch2, AliRsnPID::EType pid2);
    void   SetMotherPDG(Int_t pdg) {fMotherPDG = pdg;}

    // pair information methods
    Bool_t IsLikeSign() {return (fCharge[0] == fCharge[1]);}
    Bool_t HasEqualTypes() {return (fType[0] == fType[1]);}

    // working routines
    Double_t ComputeWeight(AliRsnDaughter *d0, AliRsnDaughter *d1);

  private:

    // pair parameters
    Int_t             fMotherPDG;  // PDG code of true mother (if known)
    Double_t          fMass[2];    // mass of particles
    Char_t            fCharge[2];  // charge of particles
    AliRsnPID::EType  fType[2];    // PID of particles

    // ROOT dictionary
    ClassDef(AliRsnPairDef, 1)
};

#endif
