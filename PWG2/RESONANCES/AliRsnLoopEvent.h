#ifndef ALIRSNLOOPEVENT_H
#define ALIRSNLOOPEVENT_H

//
// Class for computations on events
//

#include "AliRsnLoop.h"

class AliRsnLoopEvent : public AliRsnLoop {
public:

   AliRsnLoopEvent(const char *name = "default");
   AliRsnLoopEvent(const AliRsnLoopEvent &copy);
   AliRsnLoopEvent& operator=(const AliRsnLoopEvent&);
   ~AliRsnLoopEvent();
   
   virtual void       Print(Option_t *opt = "") const;
   virtual Bool_t     Init(const char *prefix, TList *list);
   virtual Int_t      DoLoop(AliRsnEvent *main, AliRsnDaughterSelector *smain = 0, AliRsnEvent *mix = 0, AliRsnDaughterSelector *smix = 0);

private:

   ClassDef(AliRsnLoopEvent, 1)
};

#endif

