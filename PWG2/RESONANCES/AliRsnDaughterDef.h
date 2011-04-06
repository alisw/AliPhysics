#ifndef ALIRSNDAUGHTERDEF_H
#define ALIRSNDAUGHTERDEF_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Single daughter definitions.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRsnDaughter.h"

class AliRsnDaughterDef : public TObject {
public:

   AliRsnDaughterDef();
   AliRsnDaughterDef(EPARTYPE type, Char_t charge = 0);
   AliRsnDaughterDef(AliRsnDaughter::ESpecies type, Char_t charge = 0);
   AliRsnDaughterDef(AliRsnDaughter::ERefType refType, Char_t charge = 0);
   AliRsnDaughterDef(const AliRsnDaughterDef &copy);
   const AliRsnDaughterDef& operator= (const AliRsnDaughterDef &copy);
   virtual ~AliRsnDaughterDef() { }

   Bool_t                   IsOnlyTrue()      const {return fOnlyTrue;}
   AliRsnDaughter::ESpecies GetPID()          const {return fPID;}
   Double_t                 GetMass()         const {return fMass;}
   Char_t                   GetChargeC()      const {return fCharge;}
   Short_t                  GetChargeS()      const {if (fCharge == '+') return 1; else if (fCharge == '-') return -1; else return 0;}
   AliRsnDaughter::ERefType GetRefType()      const {return fRefType;}
   virtual const char*      GetName()         const {return Form("%s%c", AliRsnDaughter::SpeciesName(fPID), fCharge);}
   Bool_t                   IsChargeDefined() const {return (fCharge == '+' || fCharge == '-' || fCharge == '0');}
   
   void SetOnlyTrue(Bool_t yn = kTRUE)            {fOnlyTrue = yn;}
   void SetPID(AliRsnDaughter::ESpecies pid)      {fPID = pid; fRefType = AliRsnDaughter::RefType(pid); fMass = AliRsnDaughter::SpeciesMass(pid);}
   void SetCharge(Char_t charge)                  {fCharge = charge;}
   void SetRefType(AliRsnDaughter::ERefType type) {fRefType = type;}
   
   Bool_t MatchesPID(AliRsnDaughter *daughter);
   Bool_t MatchesCharge(AliRsnDaughter *daughter);
   Bool_t MatchesRefType(AliRsnDaughter *daughter);
   Bool_t MatchesDaughter(AliRsnDaughter *daughter);
   Bool_t MatchesPDG(Int_t pdgCode)       {return (AliRsnDaughter::SpeciesPDG(fPID) == pdgCode);}
   Bool_t MatchesChargeS(Short_t charge)  {return (GetChargeS() == charge);}
   Bool_t MatchesChargeC(Char_t charge)   {return (GetChargeC() == charge);}

private:

   Bool_t                    fOnlyTrue; // fag to activate comparison of PID species
   AliRsnDaughter::ESpecies  fPID;      // PID of particles
   Double_t                  fMass;     // mass of particles (subordinate to fPID)
   Char_t                    fCharge;   // charge of particles
   AliRsnDaughter::ERefType  fRefType;  // object reference type (track/V0/cascade)

   // ROOT dictionary
   ClassDef(AliRsnDaughterDef, 1)
};

//__________________________________________________________________________________________________
inline Bool_t AliRsnDaughterDef::MatchesPID(AliRsnDaughter *daughter)
{
//
// Checks if the daughter true PID (taken from MC) matches
// that expected by this object.
// Works only if an MC is present and the fPID data member
// is set to something different from 'kUnknown'.
// If above conditions are not satisfied, it returns always kTRUE.
//

   if (fPID == AliRsnDaughter::kUnknown || daughter->GetRefMC() == 0x0) 
      return kTRUE;
   else
      return (AliRsnDaughter::SpeciesPDG(fPID) == daughter->GetPDG());
}
 
//__________________________________________________________________________________________________
inline Bool_t AliRsnDaughterDef::MatchesCharge(AliRsnDaughter *daughter)
{
//
// Checks that the daughter charge matches that expected by this object.
// Works only if the fCharge data member is set to '+', '-' or '0', 
// otherwise it accepts everything.
//

   switch (fCharge) {
      case '+': return daughter->IsPos();
      case '-': return daughter->IsNeg();
      case '0': return daughter->IsNeutral();
      default : return kTRUE;
   }
}
   
//__________________________________________________________________________________________________
inline Bool_t AliRsnDaughterDef::MatchesRefType(AliRsnDaughter *daughter)
{
//
// Checks that the daughter object type matches that expected by this object.
// Works only if the fRefType data member is different from AliRsnDaughter::kNoType,
// otherwise it accepts everything.
//

   AliRsnDaughter::ERefType type = daughter->RefType();
   
   if (fRefType != AliRsnDaughter::kNoType) 
      return (type == fRefType);
   else
      return kTRUE;
}

#endif
