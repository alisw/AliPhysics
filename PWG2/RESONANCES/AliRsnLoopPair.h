#ifndef ALIRSNLOOPPAIR_H
#define ALIRSNLOOPPAIR_H

//
// Class to combine pairs of daughters.
// See general definition inf mother class AliRsnLoop
// This is the only implementation which takes into account mixing objects.
//

#include "AliRsnLoop.h"

class TList;
class AliRsnEvent;
class AliRsnPairDef;

class AliRsnLoopPair : public AliRsnLoop {
public:

   AliRsnLoopPair(const char *name = "default", AliRsnPairDef *def = 0, Bool_t isMixed = kFALSE);
   AliRsnLoopPair(const AliRsnLoopPair &copy);
   AliRsnLoopPair& operator=(const AliRsnLoopPair &copy);
   ~AliRsnLoopPair();

   // getters
   Bool_t         IsTrueMC()     const {return fTrueMC;}
   Bool_t         IsOnlyTrue()   const {return fOnlyTrue;}
   Bool_t         IsCheckDecay() const {return fCheckDecay;}
   AliRsnPairDef* GetPairDef()         {return fPairDef;}
   AliRsnCutSet*  GetPairCuts()        {return fPairCuts;}
   AliRsnMother*  GetMother()          {return &fMother;}

   // setters (not for all members)
   void           SetPairCuts(AliRsnCutSet *cuts)      {fPairCuts = cuts;}
   void           SetTrueMC(Bool_t yn = kTRUE)         {fTrueMC = yn;}
   void           SetOnlyTrue(Bool_t onlyTrue = kTRUE) {fOnlyTrue = onlyTrue;}
   void           SetMCRefInfo(Bool_t b = kTRUE)       {fUseMCRef = b;}   
   void           SetCheckDecay(Bool_t check = kTRUE)  {fCheckDecay = check;}
   void           SetListID(Int_t i, Int_t val)        {if (i==0||i==1) fListID[i] = val;}
   void           SetRangeY(Double_t range)            {fRangeY = range;}

   // methods
   Bool_t         IsTrueMother();
   virtual void   Print(Option_t *opt = "") const;
   virtual Bool_t Init(const char *prefix, TList *list);
   virtual Int_t  DoLoop(AliRsnEvent *main, AliRsnDaughterSelector *smain, AliRsnEvent *mix = 0, AliRsnDaughterSelector *smix = 0);

protected:

   Bool_t AssignMotherAndDaughters   (AliRsnEvent *event, Int_t ipart);
   Bool_t AssignMotherAndDaughtersESD(AliRsnEvent *event, Int_t ipart);
   Bool_t AssignMotherAndDaughtersAOD(AliRsnEvent *event, Int_t ipart);
   Int_t  LoopTrueMC(AliRsnEvent *event);

   Bool_t            fTrueMC;          //  if this flag is TRUE, scan the MC for all true resonances from MC
   Bool_t            fOnlyTrue;        //  select true pairs only?
   Bool_t            fUseMCRef;        //  uses MC ref instead of REC
   Bool_t            fCheckDecay;      //  is the decay channel correct in a true pair?
   Int_t             fListID[2];       //  indexes of the two entry lists to be used
   Double_t          fRangeY;          //  range in rapidity (serves always)

   AliRsnPairDef    *fPairDef;         //  pair definition
   AliRsnCutSet     *fPairCuts;        //  collection of all cuts
   AliRsnMother      fMother;          //! mother candidate (to avoid creating it continuously)
   AliRsnDaughter    fDaughter[2];     //! daughter candidates

private:

   ClassDef(AliRsnLoopPair, 4)
};

#endif

