#ifndef ALIRSNMINIRESONANCEFINDER_H
#define ALIRSNMINIRESONANCEFINDER_H

//
// Class to identify resonance candidate pairs, which are then
// added back into the mini event to serve as daughters for
// another particle.
//

#include "TNamed.h"

#include "AliRsnMiniEvent.h"
#include "AliRsnCutSet.h"
#include "AliRsnMiniOutput.h"

//class TList;
class AliRsnMiniEvent;
//class AliRsnMiniAnalysisTask;
//class AliRsnListOutput;
//class AliRsnDaughterSelector;

class AliRsnMiniResonanceFinder : public TNamed {
public:

   AliRsnMiniResonanceFinder(const char *name = "default");
   AliRsnMiniResonanceFinder(const AliRsnMiniResonanceFinder &copy);
   AliRsnMiniResonanceFinder &operator=(const AliRsnMiniResonanceFinder &copy);
   ~AliRsnMiniResonanceFinder();

   Int_t           GetCutID(Int_t i)    const {if (i <= 0) return fCutID [0]; else return fCutID [1];}
   RSNPID          GetDaughter(Int_t i) const {if (i <= 0) return fDaughter[0]; else return fDaughter[1];}
   Double_t        GetMass(Int_t i)     const {return AliRsnDaughter::SpeciesMass(GetDaughter(i));}
   Long_t          GetPDG(Int_t i)      const {return AliRsnDaughter::SpeciesPDG(GetDaughter(i));}
   Char_t          GetCharge(Int_t i)   const {if (i <= 0) return fCharge[0]; else return fCharge[1];}
   Int_t           GetResonanceCutID()  const {return fCutIDrsn;}
   Double_t        GetResonanceMass()   const {return fResMass;}
   Long_t          GetResonancePDG()    const {return fResPDG;}
   Short_t         GetPairMode()        const {return fPairMode;}

   void            SetCutID(Int_t i, Int_t   value)   {if (i <= 0) fCutID [0] = value; else fCutID [1] = value;}
   void            SetDaughter(Int_t i, RSNPID value) {if (i <= 0) fDaughter[0] = value; else fDaughter[1] = value;}
   void            SetCharge(Int_t i, Char_t  value)  {if (i <= 0) fCharge[0] = value; else fCharge[1] = value;}
   void            SetResonanceCutID(Int_t value)     {fCutIDrsn = value;}
   void            IncrementResonanceCutID()          {fCutIDrsn++;}
   void            SetResonanceMass(Double_t mass)    {fResMass = mass;}
   void            SetResonancePDG(Long_t pdg)        {fResPDG = pdg;}
   void            SetPairCuts(AliRsnCutSet *set)     {fPairCuts = set;}
   void            SetPairMode(Short_t s)             {fPairMode = s;}

   //Int_t           ConnectTo(AliRsnMiniAnalysisTask* task);
   // ConnectTo() returns the cut ID for the resonance (an integer) which should be plugged into an AliRsnMiniOutput object.
   // However, if track cuts are subsequently added to the AliRsnMiniAnalysisTask, the cut ID for the resonance will change.
   // In any case, the function GetResonanceCutID() can be used to get the correct cut ID for the resonance.

   Int_t           RunResonanceFinder(AliRsnMiniEvent* event);
   void            FillMother(AliMCEvent* event, AliRsnMiniEvent* mini);// ESD
   void            FillMother(TClonesArray* event, AliRsnMiniEvent* mini);// AOD

private:
   
   Int_t            fCutID[2];         // ID of cut set used to select tracks
   Int_t            fCutIDrsn;         // ID of cut set that identifies the resonance candidates
   RSNPID           fDaughter[2];      // species of daughters
   Char_t           fCharge[2];        // required track charge
   Double_t         fResMass;          // nominal resonance mass
   Long_t           fResPDG;           // resonance PDG code
   AliRsnCutSet    *fPairCuts;         // cuts on the pair
   AliRsnMiniPair   fPair;             // minipair for computations
   TArrayI          fSelRes;           // list of selected resonance candidates
   TArrayI          fSel1;             // list of selected particles for definition 1
   TArrayI          fSel2;             // list of selected particles for definition 2
   Short_t          fPairMode;         // 0: use all pairs; 1: use only true pairs (MC); 2: use only false pairs (MC)

   AliRsnMiniEvent *fEvent;            // mini event to be used

   ClassDef(AliRsnMiniResonanceFinder, 1)
};

#endif
