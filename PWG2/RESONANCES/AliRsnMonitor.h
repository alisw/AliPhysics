//
// *** Class AliRsnMonitor ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef ALIRSNMonitor_H
#define ALIRSNMonitor_H

#include "TNamed.h"

#include "AliRsnDaughter.h"
#include "AliRsnDaughterDef.h"
#include "AliRsnCutSet.h"

class AliRsnMonitor : public TNamed {
public:

   AliRsnMonitor(const char *name = "default", AliRsnDaughterDef *def = 0);
   AliRsnMonitor(const AliRsnMonitor &copy);
   AliRsnMonitor& operator=(const AliRsnMonitor&);
   ~AliRsnMonitor();

   void    SetOnlyTrue(Bool_t onlyTrue = kTRUE) {fOnlyTrue = onlyTrue;}
   void    Print(Option_t *option = "") const;

   AliRsnCutSet*      GetCuts() {return &fCuts;}
   AliRsnDaughter*    GetDaughter() {return &fDaughter;}
   AliRsnDaughterDef* GetDaughterDef() {return fDaughterDef;}
   Bool_t             Fill(AliRsnDaughter *d);
   Int_t              GetCount() const {return fCount;}
   void               ResetCount() {fCount = 0;}

   virtual void       Compute();
   virtual void       Init(const char *prefix, TList *list);

protected:

   Bool_t             fOnlyTrue;        //  select true Monitors only?
   Int_t              fCount;           //  counter incremented for each added Monitor

   AliRsnDaughterDef *fDaughterDef;     //  Monitor definition (particles, charges)
   AliRsnCutSet       fCuts;            //  collection of all cuts
   AliRsnDaughter     fDaughter;        //  mother candidate (to avoid creating it continuously)

private:

   ClassDef(AliRsnMonitor, 2)
};

#endif

