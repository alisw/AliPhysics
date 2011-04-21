#ifndef ALIRSNLOOPDAUGHTER_H
#define ALIRSNLOOPDAUGHTER_H

//
// Class for computations on single daughters
//

#include "AliRsnDaughter.h"
#include "AliRsnLoop.h"

class AliRsnDaughterDef;

class AliRsnLoopDaughter : public AliRsnLoop {
public:

   AliRsnLoopDaughter(const char *name = "default", Int_t listID = 0, AliRsnDaughterDef *def = 0);
   AliRsnLoopDaughter(const AliRsnLoopDaughter &copy);
   AliRsnLoopDaughter& operator=(const AliRsnLoopDaughter&);
   ~AliRsnLoopDaughter();
   
   Int_t              GetListID()                    {return  fListID;}
   AliRsnDaughterDef* GetDef()                       {return  fDef;}
   AliRsnDaughter*    GetDaughter()                  {return &fDaughter;}
                                                     
   void               SetOnlyTrue(Bool_t yn = kTRUE) {fOnlyTrue = yn;}
   void               SetListID(Int_t i)             {fListID = i;}
   void               SetDef(AliRsnDaughterDef *def) {fDef = def;}
   
   virtual void       Print(Option_t *opt = "") const;
   virtual Bool_t     Init(const char *prefix, TList *list);
   virtual Int_t      DoLoop(AliRsnEvent *main, AliRsnDaughterSelector *smain, AliRsnEvent *mix = 0, AliRsnDaughterSelector *smix = 0);

protected:

   Bool_t             fOnlyTrue;   //  for selecting only true particles
   Int_t              fListID;     //  index of entry list to use
   AliRsnDaughterDef *fDef;        //  definition for selection
   AliRsnDaughter     fDaughter;   //! daughter temporary member

private:

   ClassDef(AliRsnLoopDaughter, 3)
};

#endif

