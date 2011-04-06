#ifndef ALIRSNLOOP_H
#define ALIRSNLOOP_H

//
// Base class object for all computation devices.
//

#include "TNamed.h"
#include "TObjArray.h"

#include "AliRsnListOutput.h"
#include "AliRsnValue.h"
#include "AliRsnCutSet.h"

class TList;
class AliRsnEvent;
class AliRsnListOutput;
class AliRsnDaughterSelector;

class AliRsnLoop : public TNamed {
public:

   enum EOut {
      kH1,
      kHSparse,
      kNtuple
   };

   AliRsnLoop(const char *name = "default", Bool_t isMixed = kFALSE);
   AliRsnLoop(const AliRsnLoop &copy);
   AliRsnLoop& operator=(const AliRsnLoop&);
   ~AliRsnLoop();
   
   void           SetMixed(Bool_t yn = kTRUE)     {fIsMixed = yn;}
   void           SetEventCuts(AliRsnCutSet *set) {fEventCuts = set;}
   Bool_t         IsMixed() const                 {return fIsMixed;}
   AliRsnCutSet*  GetEventCuts()                  {return fEventCuts;}
   Bool_t         OkEvent(AliRsnEvent *rsn);
   
   virtual void   AddOutput(TObject *output);
   virtual void   Print(Option_t *option = "") const;
   virtual Bool_t Init(const char *prefix, TList *list);
   virtual Int_t  DoLoop(AliRsnEvent *main, AliRsnDaughterSelector *smain, AliRsnEvent *mix = 0, AliRsnDaughterSelector *smix = 0);

protected:

   Bool_t         fIsMixed;    //  flag to know if the loop works with event mixing
   AliRsnCutSet  *fEventCuts;  //  event cuts
   TClonesArray   fOutputs;    //  output object definitions

private:

   ClassDef(AliRsnLoop, 1)
};

#endif
