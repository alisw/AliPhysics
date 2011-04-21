#ifndef ALIRSNDAUGHTERDEF_H
#define ALIRSNDAUGHTERDEF_H

//
//  Definition for a single candidate daughter.
//  They can be chosen among the following possibilities:
//
//  -- particle species, chosen in the enum AliRsnDaughter::ESpecies;
//     this option does two things:
//     -- when possible (needs MC) and required, allow to check is a daughter
//        is of the defined species (e.g.: for selecting true daughters of a resonance)
//     -- defines the mass to be assigned to this object, for any purpose
//        (e.g.: to compute its 4-momentum to be used for a resonance mass)
//
//  -- track charge, which can be '+', '-', '0', 
//     -- any other char leaves this undefined (use only in daughter monitor loops)
//     -- when doing resonance analysis, or when RSN Input handler needs to be used,
//        this must always be defined
//     
//  -- object type (track/V0/cascade)
//     -- could be needed to select tracks when particle species is not specified
//     -- works only in single daughter loops
//     

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

   AliRsnDaughter::ESpecies GetPID()          const {return fPID;}
   Double_t                 GetMass()         const {return fMass;}
   Char_t                   GetChargeC()      const {return fCharge;}
   Short_t                  GetChargeS()      const {if (fCharge == '+') return 1; else if (fCharge == '-') return -1; else return 0;}
   AliRsnDaughter::ERefType GetRefType()      const {return fRefType;}
   virtual const char*      GetName()         const {return Form("%s%c", AliRsnDaughter::SpeciesName(fPID), fCharge);}
   Bool_t                   IsChargeDefined() const {return (fCharge == '+' || fCharge == '-' || fCharge == '0');}
   
   void SetPID(AliRsnDaughter::ESpecies pid)      {fPID = pid; fRefType = AliRsnDaughter::RefType(pid); fMass = AliRsnDaughter::SpeciesMass(pid);}
   void SetCharge(Char_t charge)                  {fCharge = charge;}
   void SetRefType(AliRsnDaughter::ERefType type) {fRefType = type;}
   
   Bool_t MatchesPID(AliRsnDaughter *daughter);
   Bool_t MatchesCharge(AliRsnDaughter *daughter);
   Bool_t MatchesRefType(AliRsnDaughter *daughter);
   Bool_t MatchesPDG(Int_t pdgCode)                  {return (AliRsnDaughter::SpeciesPDG(fPID) == pdgCode);}
   Bool_t MatchesChargeS(Short_t charge)             {return (GetChargeS() == charge);}
   Bool_t MatchesChargeC(Char_t charge)              {return (GetChargeC() == charge);}

private:

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
// Checks if the passed daughter true particle type 
// matches the species defined in fPID data member,
// by comparing the corresponding PDG codes of both in absolute value.
// Returns kTRUE when the two codes match, and kFALSE when:
// -- PDG codes don't match
// -- fPID is not set to a well-defined species (AliRsnDaughter::kUnknown)
// -- passed daughter has not MC information
//
   if (fPID == AliRsnDaughter::kUnknown) {
      AliError("This DaughterDef has undefined species: cannot check PDG matching with passed arg");
      return kFALSE;
   } else if (!daughter->GetRefMC()) {
      AliError("The passed argument has NULL MC pointer: cannot check PDG matching with this DaughterDef");
      return kFALSE;
   }
   
   return MatchesPDG(daughter->GetPDGAbs());
}
 
//__________________________________________________________________________________________________
inline Bool_t AliRsnDaughterDef::MatchesCharge(AliRsnDaughter *daughter)
{
//
// Checks if the passed daughter charge matches that defined in fCharge data member.
// If fCharge is not initialized to '+', '-' or '0', that is interpreted as if the user
// does not want to select in charge, and then the function returns always kTRUE.
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
// Checks that the daughter object type matches that defined in fRefType data member.
// If fRefType is initialized to AliRsnDaughter::kNoType, that is interpreted as if the user
// does not want to select in type and then the function returns always kTRUE.
//

   AliRsnDaughter::ERefType type = daughter->RefType();
   
   if (fRefType != AliRsnDaughter::kNoType) 
      return (type == fRefType);
   else
      return kTRUE;
}

#endif
