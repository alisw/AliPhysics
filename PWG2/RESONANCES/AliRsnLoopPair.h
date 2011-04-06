#ifndef ALIRSNLOOPPAIR_H
#define ALIRSNLOOPPAIR_H

//
// Class to combine pairs of daughters.
//

#include "AliRsnPairDef.h"
#include "AliRsnLoop.h"

class TList;

class AliRsnLoopPair : public AliRsnLoop {
public:

   AliRsnLoopPair(const char *name = "default", AliRsnPairDef *def = 0, Bool_t isMixed = kFALSE);
   AliRsnLoopPair(const AliRsnLoopPair &copy);
   AliRsnLoopPair& operator=(const AliRsnLoopPair&);
   ~AliRsnLoopPair();

   // getters
   Bool_t         IsOnlyTrue()   const {return fOnlyTrue;}
   Bool_t         IsCheckDecay() const {return fCheckDecay;}
   AliRsnPairDef* GetPairDef()         {return fPairDef;}
   AliRsnCutSet*  GetPairCuts()        {return fPairCuts;}
   AliRsnMother*  GetMother()          {return &fMother;}

   // setters (not for all members)
   void           SetPairCuts(AliRsnCutSet *cuts)      {fPairCuts = cuts;}
   void           SetOnlyTrue(Bool_t onlyTrue = kTRUE) {fOnlyTrue = onlyTrue;}
   void           SetCheckDecay(Bool_t check = kTRUE)  {fCheckDecay = check;}
   void           SetListID(Int_t i, Int_t val)        {if (i==0||i==1) fListID[i] = val;}

   // methods
   Bool_t         MotherOK();
   virtual void   Print(Option_t *opt = "") const;
   virtual Bool_t Init(const char *prefix, TList *list);
   virtual Int_t  DoLoop(AliRsnEvent *main, AliRsnDaughterSelector *smain, AliRsnEvent *mix = 0, AliRsnDaughterSelector *smix = 0);

protected:

   Bool_t            fOnlyTrue;        //  select true pairs only?
   Bool_t            fCheckDecay;      //  is the decay channel correct in a true pair?
   Int_t             fListID[2];       //  indexes of the two entry lists to be used

   AliRsnPairDef    *fPairDef;         //  pair definition
   AliRsnCutSet     *fPairCuts;        //  collection of all cuts
   AliRsnMother      fMother;          //! mother candidate (to avoid creating it continuously)
   AliRsnDaughter    fDaughter[2];     //! daughter candidates

private:

   ClassDef(AliRsnLoopPair, 2)
};

#endif

