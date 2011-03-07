//
// *** Class AliRsnPair ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef ALIRSNPAIR_H
#define ALIRSNPAIR_H

#include "TNamed.h"

#include "AliRsnPairDef.h"
#include "AliRsnCutManager.h"
#include "AliRsnMother.h"

class TList;

class AliRsnPair : public TNamed {
public:

   AliRsnPair(const char *name = "default", AliRsnPairDef *def = 0);
   AliRsnPair(const AliRsnPair &copy);
   AliRsnPair& operator=(const AliRsnPair&);
   ~AliRsnPair();
   
   // getters
   Bool_t            IsOnlyTrue()   const {return fOnlyTrue;}
   Bool_t            IsCheckDecay() const {return fCheckDecay;}
   Bool_t            IsMixed()      const {return fIsMixed;}
   Int_t             GetCount()     const {return fCount;}
   AliRsnPairDef*    GetPairDef()         {return fPairDef;}
   AliRsnCutManager* GetCutManager()      {return &fCutManager;}
   AliRsnMother*     GetMother()          {return &fMother;}
   
   // shortcuts to data-member getters
   AliRsnCutSet*     GetCommonDaughterCuts() {return fCutManager.GetCommonDaughterCuts();}
   AliRsnCutSet*     GetDaughter1Cuts()      {return fCutManager.GetDaughter1Cuts();}
   AliRsnCutSet*     GetDaughter2Cuts()      {return fCutManager.GetDaughter2Cuts();}
   AliRsnCutSet*     GetMotherCuts()         {return fCutManager.GetMotherCuts();}
                     
   // setters (not for all members)
   void              SetOnlyTrue(Bool_t onlyTrue = kTRUE) {fOnlyTrue = onlyTrue;}
   void              SetCheckDecay(Bool_t check = kTRUE)  {fCheckDecay = check;}
   void              SetMixed(Bool_t doit = kTRUE)        {fIsMixed = doit;}
   void              SetCount(Int_t count)                {fCount = count;}
   void              ResetCount()                         {fCount = 0;}

   // methods
   Bool_t            Fill(AliRsnDaughter *d0, AliRsnDaughter *d1, Bool_t refFirst = kTRUE);
   virtual void      Print(Option_t *option = "") const;
   virtual void      Compute();
   virtual void      Init(const char *prefix, TList *list);

protected:

   Bool_t            fOnlyTrue;        //  select true pairs only?
   Bool_t            fCheckDecay;      //  is the decay channel correct in a true pair?
   Bool_t            fIsMixed;         //  is this an event-mixing?
   Int_t             fCount;           //  counter incremented for each added pair

   AliRsnPairDef    *fPairDef;         //  pair definition (particles, charges)
   AliRsnCutManager  fCutManager;      //  collection of all cuts
   AliRsnMother      fMother;          //  mother candidate (to avoid creating it continuously)

private:

   ClassDef(AliRsnPair, 2)
};

#endif

