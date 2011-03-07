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
   
   // getters
   Bool_t             IsOnlyTrue()   const {return fOnlyTrue;}
   Int_t              GetCount()     const {return fCount;}
   AliRsnDaughterDef* GetDaughterDef()     {return fDaughterDef;}
   AliRsnCutSet*      GetCuts()            {return &fCuts;}
   AliRsnDaughter*    GetDaughter()        {return fDaughter;}
   
   // setters (not for all members)
   void               SetOnlyTrue(Bool_t onlyTrue = kTRUE) {fOnlyTrue = onlyTrue;}
   void               SetCount(Int_t count)                {fCount = count;}
   void               ResetCount()                         {fCount = 0;}
   
   // methods
   Bool_t             Fill(AliRsnDaughter *d);
   virtual void       Print(Option_t *option = "") const;
   virtual void       Compute();
   virtual void       Init(const char *prefix, TList *list);

protected:

   Bool_t             fOnlyTrue;        //  select true Monitors only?
   Int_t              fCount;           //  counter incremented for each added Monitor

   AliRsnDaughterDef *fDaughterDef;     //  Monitor definition (particles, charges)
   AliRsnCutSet       fCuts;            //  collection of all cuts
   AliRsnDaughter    *fDaughter;        //! pointer to daughter

private:

   ClassDef(AliRsnMonitor, 2)
};

#endif

