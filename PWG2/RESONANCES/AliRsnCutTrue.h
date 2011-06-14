#ifndef ALIRSNCUTTRUE_H
#define ALIRSNCUTTRUE_H

//
// This cut selects the AliRsnDaughter objects pointing
// to tracks with a well defined true particle species,
// defined through its PDG code or species according the
// enumeration defined in AliRsnDaughter class.
// ---
// Using this cut on data results in no tracks passing it.
//

#include "AliRsnDaughter.h"
#include "AliRsnCut.h"

class AliRsnCutTrue : public AliRsnCut {
   
public:

   AliRsnCutTrue(const char *name, Int_t pdg);
   AliRsnCutTrue(const char *name, AliRsnDaughter::ESpecies species);
   AliRsnCutTrue(const AliRsnCutTrue& copy);
   AliRsnCutTrue& operator=(const AliRsnCutTrue &copy);
   virtual ~AliRsnCutTrue() { }
   
   virtual Bool_t IsSelected(TObject *obj);
   
private:

   ClassDef(AliRsnCutTrue,1)    // AliRsnCutTrue class

};

#endif
