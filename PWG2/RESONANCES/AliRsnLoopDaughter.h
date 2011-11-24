#ifndef ALIRSNLOOPDAUGHTER_H
#define ALIRSNLOOPDAUGHTER_H

//
// Computator for single daughters.
// Implements a simple loop on tracks from one of the entry lists
// filled by the task AliRsnInputHandler, adding a check on their
// definition specified in the daughter def.
//

#include "AliRsnDaughter.h"
#include "AliRsnLoop.h"

class AliRsnDaughterDef;

class AliRsnLoopDaughter : public AliRsnLoop {
public:

   AliRsnLoopDaughter(const char *name = "default", Int_t listID = 0, AliRsnDaughterDef *def = 0);
   AliRsnLoopDaughter(const AliRsnLoopDaughter &copy);
   AliRsnLoopDaughter& operator=(const AliRsnLoopDaughter &copy);
   ~AliRsnLoopDaughter();
   
   Int_t              GetListID() const              {return  fListID;}
   AliRsnDaughterDef* GetDef()                       {return  fDef;}
   AliRsnDaughter*    GetDaughter()                  {return &fDaughter;}

   void               SetTrueMC(Bool_t yn = kTRUE)   {fTrueMC = yn;}
   void               SetOnlyTrue(Bool_t yn = kTRUE) {fOnlyTrue = yn;}
   void               SetMCRefInfo(Bool_t b = kTRUE) {fUseMCRef = b;}  
   void               SetListID(Int_t i)             {fListID = i;}
   void               SetDef(AliRsnDaughterDef *def) {fDef = def;}
   
   virtual void       Print(Option_t *opt = "") const;
   virtual Bool_t     Init(const char *prefix, TList *list);
   virtual Int_t      DoLoop(AliRsnEvent *main, AliRsnDaughterSelector *smain, AliRsnEvent *mix = 0, AliRsnDaughterSelector *smix = 0);

protected:

   Int_t LoopTrueMC(AliRsnEvent *rsn);

   Bool_t             fTrueMC;     //  if this flag is TRUE, scan the MC for all true resonances from MC
   Bool_t             fOnlyTrue;   //  for selecting only true particles
   Bool_t             fUseMCRef;   //  uses MC ref instead of REC
   Int_t              fListID;     //  index of entry list to use
   AliRsnDaughterDef *fDef;        //  definition for selection
   AliRsnDaughter     fDaughter;   //! daughter temporary member

private:

   ClassDef(AliRsnLoopDaughter, 4)
};

#endif

