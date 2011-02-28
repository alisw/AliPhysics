#ifndef ALIRSNPAIRDEF_H
#define ALIRSNPAIRDEF_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Resonance decay tree definition.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRsnDaughter.h"
#include "AliRsnDaughterDef.h"

class AliRsnPairDef : public TObject {
public:

   AliRsnPairDef();
   AliRsnPairDef(EPARTYPE type1, Char_t sign1, EPARTYPE type2, Char_t sign2, Int_t motherPDG = 0, Double_t motherMass = 0.0);
   AliRsnPairDef(AliRsnDaughter::ESpecies type1, Char_t sign1, AliRsnDaughter::ESpecies type2, Char_t sign2, Int_t motherPDG = 0, Double_t motherMass = 0.0);
   AliRsnPairDef(const AliRsnPairDef &copy);
   const AliRsnPairDef& operator= (const AliRsnPairDef &copy);
   virtual ~AliRsnPairDef() { }

   // getters
   Int_t                    GetMotherPDG()  const {return fMotherPDG;}
   Double_t                 GetMotherMass() const {return fMotherMass;}
   const char*              GetPairName()   const;
   virtual const char*	    GetName()       const {return GetPairName();}
   
   AliRsnDaughterDef*       GetDef1()           {return &fDef1;}
   Double_t                 GetMass1()    const {return fDef1.GetMass();}
   Char_t                   GetChargeC1() const {return fDef1.GetChargeC();}
   Short_t                  GetChargeS1() const {return fDef1.GetChargeS();}
   AliRsnDaughter::ESpecies GetPID1()     const {return fDef1.GetPID();}
   AliRsnDaughter::ERefType GetRefType1() const {return fDef1.GetRefType();}
   
   AliRsnDaughterDef*       GetDef2()           {return &fDef2;}
   Double_t                 GetMass2()    const {return fDef2.GetMass();}
   Char_t                   GetChargeC2() const {return fDef2.GetChargeC();}
   Short_t                  GetChargeS2() const {return fDef2.GetChargeS();}
   AliRsnDaughter::ESpecies GetPID2()     const {return fDef2.GetPID();}
   AliRsnDaughter::ERefType GetRefType2() const {return fDef2.GetRefType();}
   
   AliRsnDaughterDef*       GetDef(Int_t i)           {if (i<1) return GetDef1(); else return GetDef2();}
   Double_t                 GetMass(Int_t i)    const {if (i<1) return GetMass1(); else return GetMass2();}
   Char_t                   GetChargeC(Int_t i) const {if (i<1) return GetChargeC1(); else return GetChargeS2();}
   Short_t                  GetChargeS(Int_t i) const {if (i<1) return GetChargeS1(); else return GetChargeS2();}
   AliRsnDaughter::ESpecies GetPID(Int_t i)     const {if (i<1) return GetPID1(); else return GetPID2();}
   AliRsnDaughter::ERefType GetRefType(Int_t i) const {if (i<1) return GetRefType1(); else return GetRefType2();}

   // setters
   void SetMotherPDG(Int_t pdg)                              {fMotherPDG = pdg;}
   void SetMotherMass(Double_t mass)                         {fMotherMass = mass;}
   void SetDef1(EPARTYPE pid, Char_t charge)                 {fDef1.SetCharge(charge); fDef1.SetPID(AliRsnDaughter::FromAliPID(pid));}
   void SetDef2(EPARTYPE pid, Char_t charge)                 {fDef2.SetCharge(charge); fDef2.SetPID(AliRsnDaughter::FromAliPID(pid));}
   void SetDef1(AliRsnDaughter::ESpecies pid, Char_t charge) {fDef1.SetCharge(charge); fDef1.SetPID(pid);}
   void SetDef2(AliRsnDaughter::ESpecies pid, Char_t charge) {fDef2.SetCharge(charge); fDef2.SetPID(pid);}
   
   void SetDefs(EPARTYPE pid1, Char_t ch1, AliPID::EParticleType pid2, Char_t ch2)                 {SetDef1(pid1, ch1); SetDef2(pid2, ch2);}
   void SetDefs(AliRsnDaughter::ESpecies pid1, Char_t ch1, AliPID::EParticleType pid2, Char_t ch2) {SetDef1(pid1, ch1); SetDef2(pid2, ch2);}

   // checkers
   Bool_t IsLikeSign()  const {return (GetChargeC1() == GetChargeC2());}
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
