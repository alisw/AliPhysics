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
#include "AliRsnDaughter.h"

class AliRsnDaughter;

class AliRsnPairDef : public TObject
{
  public:

    AliRsnPairDef();
    AliRsnPairDef(AliPID::EParticleType type1, Char_t sign1, AliPID::EParticleType type2, Char_t sign2, Int_t motherPDG = 0, Double_t motherMass = 0.0);
    AliRsnPairDef(const AliRsnPairDef &copy);
    const AliRsnPairDef& operator= (const AliRsnPairDef &copy);
    virtual ~AliRsnPairDef() { }

    // getters
    Int_t                    GetMotherPDG()  const {return fMotherPDG;}
    Double_t                 GetMotherMass() const {return fMotherMass;}
    Double_t                 GetMass(Int_t i) const {if (i>=0&&i<2) return fMass[i]; else return 0.0;}
    Char_t                   GetCharge(Int_t i) const {if (i>=0&&i<2) return fCharge[i]; else return 0;}
    Short_t                  GetChargeShort(Int_t i) const {if (GetCharge(i) == '+') return 1; else if (GetCharge(i) == '-') return -1; else return 0;}
    AliPID::EParticleType    GetPID(Int_t i) const {if (i>=0&&i<2) return fPID[i]; else return AliPID::kUnknown;}
    AliRsnDaughter::ERefType GetDaughterType(Int_t i) {if (i>=0&&i<2) return fDaughterType[i]; else return AliRsnDaughter::kNoType;}
    const char*              GetPairName() const;

    // setters
    void   SetMotherPDG(Int_t pdg) {fMotherPDG = pdg;}
    void   SetMotherMass(Double_t mass) {fMotherMass = mass;}
    Bool_t SetDaughter(Int_t i, AliPID::EParticleType pid, Char_t charge = 0);
    Bool_t SetDaughters(AliPID::EParticleType type1, Char_t sign1, AliPID::EParticleType type2, Char_t sign2);

    // pair information methods
    Bool_t IsLikeSign() const {return (fCharge[0] == fCharge[1]);}
    Bool_t HasEqualPID() const {return (fPID[0] == fPID[1]);}

  private:

    // pair parameters
    Double_t                  fMotherMass;      // nominal mass of true mother
    Int_t                     fMotherPDG;       // PDG code of true mother (if known)
    
    Double_t                  fMass[2];         // mass of particles
    Char_t                    fCharge[2];       // charge of particles
    AliPID::EParticleType     fPID[2];          // PID of particles
    AliRsnDaughter::ERefType  fDaughterType[2]; // object type (track/V0)

    // ROOT dictionary
    ClassDef(AliRsnPairDef, 1)
};

#endif
