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

#include "AliPID.h"

class AliRsnDaughter;

class AliRsnPairDef : public TObject
{
  public:

    AliRsnPairDef();
    AliRsnPairDef(Char_t sign1, AliPID::EParticleType type1,
                  Char_t sign2, AliPID::EParticleType type2, Int_t motherPDG = 0);
    AliRsnPairDef(AliPID::EParticleType type1, Char_t sign1,
                  AliPID::EParticleType type2, Char_t sign2, Int_t motherPDG = 0);
    AliRsnPairDef(const AliRsnPairDef &copy);
    const AliRsnPairDef& operator= (const AliRsnPairDef &copy);
    virtual ~AliRsnPairDef() { }

    // getters
    Char_t                GetCharge(Int_t i) const {if (i>=0&&i<2) return fCharge[i]; else return 0;}
    Short_t               GetChargeShort(Int_t i) const {if (GetCharge(i) == '+') return 1; else if (GetCharge(i) == '-') return -1; else return 0;}
    AliPID::EParticleType GetType(Int_t i) const {if (i>=0&&i<2) return fType[i]; else return AliPID::kUnknown;}
    Double_t              GetMass(Int_t i) const {if (i>=0&&i<2) return fMass[i]; else return 0.0;}
    Int_t                 GetMotherPDG() const {return fMotherPDG;}
    TString               GetPairName() const;

    // setters
    Bool_t SetPairElement(Int_t i, Char_t charge, AliPID::EParticleType pid);
    Bool_t SetPair(Char_t ch1, AliPID::EParticleType pid1, Char_t ch2, AliPID::EParticleType pid2);
    void   SetMotherPDG(Int_t pdg) {fMotherPDG = pdg;}

    // pair information methods
    Bool_t IsLikeSign() const {return (fCharge[0] == fCharge[1]);}
    Bool_t HasEqualTypes() const {return (fType[0] == fType[1]);}

  private:

    // pair parameters
    Int_t                  fMotherPDG;  // PDG code of true mother (if known)
    Double_t               fMass[2];    // mass of particles
    Char_t                 fCharge[2];  // charge of particles
    AliPID::EParticleType  fType[2];    // PID of particles

    // ROOT dictionary
    ClassDef(AliRsnPairDef, 1)
};

#endif
