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
   AliRsnPairDef(EPARTYPE type1, Char_t ch1, EPARTYPE type2, Char_t ch2, Int_t motherPDG = 0, Double_t motherMass = 0.0);
   AliRsnPairDef(AliRsnDaughter::ESpecies type1, Char_t ch1, AliRsnDaughter::ESpecies type2, Char_t ch2, Int_t motherPDG = 0, Double_t motherMass = 0.0);
   AliRsnPairDef(const AliRsnPairDef &copy);
   const AliRsnPairDef& operator= (const AliRsnPairDef &copy);
   virtual ~AliRsnPairDef() { }

   virtual const char*	    GetName()       const {return Form("%s_%s", fDef1.GetName(), fDef2.GetName());}
   Int_t                    GetMotherPDG()  const {return fMotherPDG;}
   Double_t                 GetMotherMass() const {return fMotherMass;}
   AliRsnDaughterDef&       GetDef1()             {return fDef1;}
   AliRsnDaughterDef&       GetDef2()             {return fDef2;}
   AliRsnDaughterDef&       GetDef(Int_t i)       {if (i<1) return GetDef1(); else return GetDef2();}

   void SetMotherPDG(Int_t pdg)                 {fMotherPDG = pdg;}
   void SetMotherMass(Double_t mass)            {fMotherMass = mass;}
   void SetDef1(AliRsnDaughterDef *def)         {if (def) fDef1 = (*def);}
   void SetDef2(AliRsnDaughterDef *def)         {if (def) fDef2 = (*def);}
   void SetDef(Int_t i, AliRsnDaughterDef *def) {if (!def) return; if (i<1) fDef1 = (*def); else fDef2 = (*def);}

   Bool_t IsLikeSign()  const {return (fDef1.GetChargeC() == fDef2.GetChargeC());}
   Bool_t HasEqualPID() const {return (fDef1.GetPID() == fDef2.GetPID());}

private:

   Double_t          fMotherMass;  // nominal mass of true mother
   Int_t             fMotherPDG;   // PDG code of true mother (if known)
   AliRsnDaughterDef fDef1;        // definitions for daughter #1 (see class)
   AliRsnDaughterDef fDef2;        // definitions for daughter #2 (see class)

   // ROOT dictionary
   ClassDef(AliRsnPairDef, 1)
};

#endif
