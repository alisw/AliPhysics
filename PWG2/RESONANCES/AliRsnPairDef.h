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
#include "AliRsnDaughterDef.h"

class AliRsnMother;

class AliRsnPairDef : public TObject {
public:

   AliRsnPairDef();
   AliRsnPairDef(EPARTYPE type1, Char_t sign1, EPARTYPE type2, Char_t sign2, Int_t motherPDG = 0, Double_t motherMass = 0.0);
   AliRsnPairDef(AliRsnDaughterDef def1, AliRsnDaughterDef def2);
   AliRsnPairDef(const AliRsnPairDef &copy);
   const AliRsnPairDef& operator= (const AliRsnPairDef &copy);
   virtual ~AliRsnPairDef() { }

   // getters
   Int_t                    GetMotherPDG()  const {return fMotherPDG;}
   Double_t                 GetMotherMass() const {return fMotherMass;}
   const char*              GetPairName()   const;
   
   AliRsnDaughterDef*       GetDef1()                {return &fDef1;}
   Double_t                 GetMass1()         const {return fDef1.GetMass();}
   Char_t                   GetCharge1()       const {return fDef1.GetCharge();}
   Short_t                  GetChargeShort1()  const {return fDef1.GetChargeShort();}
   AliPID::EParticleType    GetPID1()          const {return fDef1.GetPID();}
   AliRsnDaughter::ERefType GetDaughterType1() const {return fDef1.GetDaughterType();}
   
   AliRsnDaughterDef*       GetDef2()                {return &fDef2;}
   Double_t                 GetMass2()         const {return fDef2.GetMass();}
   Char_t                   GetCharge2()       const {return fDef2.GetCharge();}
   Short_t                  GetChargeShort2()  const {return fDef2.GetChargeShort();}
   AliPID::EParticleType    GetPID2()          const {return fDef2.GetPID();}
   AliRsnDaughter::ERefType GetDaughterType2() const {return fDef2.GetDaughterType();}
   
   AliRsnDaughterDef*       GetDef(Int_t i)                {if (i<1) return GetDef1(); else return GetDef2();}
   Double_t                 GetMass(Int_t i)         const {if (i<1) return GetMass1(); else return GetMass2();}
   Char_t                   GetCharge(Int_t i)       const {if (i<1) return GetCharge1(); else return GetCharge2();}
   Short_t                  GetChargeShort(Int_t i)  const {if (i<1) return GetChargeShort1(); else return GetChargeShort2();}
   AliPID::EParticleType    GetPID(Int_t i)          const {if (i<1) return GetPID1(); else return GetPID2();}
   AliRsnDaughter::ERefType GetDaughterType(Int_t i) const {if (i<1) return GetDaughterType1(); else return GetDaughterType2();}

   // setters
   void   SetMotherPDG(Int_t pdg)                           {fMotherPDG = pdg;}
   void   SetMotherMass(Double_t mass)                      {fMotherMass = mass;}
   Bool_t SetDef1(AliPID::EParticleType pid, Char_t charge) {return fDef1.SetDaughter(pid, charge);}
   Bool_t SetDef2(AliPID::EParticleType pid, Char_t charge) {return fDef2.SetDaughter(pid, charge);}
   Bool_t SetDefs(AliPID::EParticleType pid1, Char_t ch1, AliPID::EParticleType pid2, Char_t ch2) {return (SetDef1(pid1, ch1) && SetDef2(pid2, ch2));}

   // checkers
   Bool_t IsLikeSign()  const {return (GetCharge1() == GetCharge2());}
   Bool_t HasEqualPID() const {return (GetPID1() == GetPID2());}

private:

   Double_t          fMotherMass;  // nominal mass of true mother
   Int_t             fMotherPDG;   // PDG code of true mother (if known)
   AliRsnDaughterDef fDef1;        // definitions for daughter #1 (see class)
   AliRsnDaughterDef fDef2;        // definitions for daughter #2 (see class)

   // ROOT dictionary
   ClassDef(AliRsnPairDef, 1)
};

#endif
