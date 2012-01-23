#ifndef ALIDIELECTRONSIGNALMC_H
#define ALIDIELECTRONSIGNALMC_H

#include <TNamed.h>


/*
   Monte Carlo signal definition:
      Leg #1  <-- Mother #1  <--  Grandmother #1
                      |
      Leg #2  <-- Mother #2  <--  Grandmother #2
   All particles can be classified as:
     1. Primary   - particle originating in the physics event
     2. Secondary - particle created during the GEANT propagation due to interaction of primaries with the material
     3. Direct    - particle directly created in the collision (has no mother)
     4. Secondary - particle which is the product of the decay or reinteraction of another particle
   The 2 legs can originate from the same or different mother particles.
*/


//__________________________________________________________________
class AliDielectronSignalMC : public TNamed {
  
 public:
  enum EBranchRelation {kUndefined=0, kSame, kDifferent};
  enum ESource {kDontCare=0, kPrimary, kSecondary, kDirect, kDecayProduct};
  
  AliDielectronSignalMC();
  AliDielectronSignalMC(const Char_t* name, const Char_t* title);
  virtual ~AliDielectronSignalMC();
  
  void SetLegPDGs(Int_t pdg1, Int_t pdg2)                 {fLeg1 = pdg1; fLeg2 = pdg2;}
  void SetMotherPDGs(Int_t pdg1, Int_t pdg2)              {fMother1 = pdg1; fMother2 = pdg2;}
  void SetGrandMotherPDGs(Int_t pdg1, Int_t pdg2)         {fGrandMother1 = pdg1; fGrandMother2 = pdg2;}
  void SetLegSources(ESource s1, ESource s2)              {fLeg1Source = s1; fLeg2Source = s2;}
  void SetMotherSources(ESource s1, ESource s2)           {fMother1Source = s1; fMother2Source = s2;}
  void SetGrandMotherSources(ESource s1, ESource s2)      {fGrandMother1Source = s1; fGrandMother2Source = s2;}
  void SetCheckBothChargesLegs(Bool_t flag1, Bool_t flag2) {fCheckBothChargesLeg1 = flag1; fCheckBothChargesLeg2 = flag2;}
  void SetCheckBothChargesMothers(Bool_t flag1, Bool_t flag2) {fCheckBothChargesMother1 = flag1; fCheckBothChargesMother2 = flag2;}
  void SetCheckBothChargesGrandMothers(Bool_t flag1, Bool_t flag2) {fCheckBothChargesGrandMother1 = flag1; fCheckBothChargesGrandMother2 = flag2;}
  void SetMothersRelation(EBranchRelation relation)       {fMothersRelation = relation;}
  void SetFillPureMCStep(Bool_t fill=kTRUE)               {fFillPureMCStep = fill;}

  Int_t GetLegPDG(Int_t branch) const                     {return (branch==1 ? fLeg1 : fLeg2);}
  Int_t GetMotherPDG(Int_t branch) const                  {return (branch==1 ? fMother1 : fMother2);}
  Int_t GetGrandMotherPDG(Int_t branch) const             {return (branch==1 ? fGrandMother1 : fGrandMother2);}
  ESource GetLegSource(Int_t branch) const                {return (branch==1 ? fLeg1Source : fLeg2Source);}
  ESource GetMotherSource(Int_t branch) const             {return (branch==1 ? fMother1Source : fMother2Source);}
  ESource GetGrandMotherSource(Int_t branch) const        {return (branch==1 ? fGrandMother1Source : fGrandMother2Source);}
  Bool_t GetCheckBothChargesLegs(Int_t branch) const      {return (branch==1 ? fCheckBothChargesLeg1 : fCheckBothChargesLeg2);}
  Bool_t GetCheckBothChargesMothers(Int_t branch) const   {return (branch==1 ? fCheckBothChargesMother1 : fCheckBothChargesMother2);}
  Bool_t GetCheckBothChargesGrandMothers(Int_t branch) const   {return (branch==1 ? fCheckBothChargesGrandMother1 : fCheckBothChargesGrandMother2);}
  EBranchRelation GetMothersRelation() const              {return fMothersRelation;}
  Bool_t GetFillPureMCStep() const                        {return fFillPureMCStep;}

 private:
  Int_t fLeg1;                        // leg 1 PDG
  Int_t fLeg2;                        // leg 2 PDG
  Int_t fMother1;                     // mother 1 PDG
  Int_t fMother2;                     // mother 2 PDG
  Int_t fGrandMother1;                // grandmother 1 PDG
  Int_t fGrandMother2;                // grandmother 2 PDG
    
  ESource fLeg1Source;                // leg 1 source
  ESource fLeg2Source;                // leg 2 source
  ESource fMother1Source;             // mother 1 source
  ESource fMother2Source;             // mother 2 source
  ESource fGrandMother1Source;        // grandmother 1 source
  ESource fGrandMother2Source;        // grandmother 2 source

  Bool_t fCheckBothChargesLeg1;         // check both charges of the legs pdg
  Bool_t fCheckBothChargesLeg2;         //                leg2
  Bool_t fCheckBothChargesMother1;      //                mother 1
  Bool_t fCheckBothChargesMother2;      //                mother 2
  Bool_t fCheckBothChargesGrandMother1; //              grand mother 1
  Bool_t fCheckBothChargesGrandMother2; //              grand mother 2
  
  EBranchRelation fMothersRelation;   // mother 1&2 relation (same, different or whatever)
  
  Bool_t fFillPureMCStep;             // check and fill the pure MC step
  
  ClassDef(AliDielectronSignalMC,1);
};

#endif
