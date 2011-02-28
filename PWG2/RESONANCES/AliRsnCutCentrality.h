#ifndef ALIRSNCUTCENTRALITY_H
#define ALIRSNCUTCENTRALITY_H

////////////////////////////////////////////////////////////////////////////////
//
//  Centrality cut
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRsnCut.h"

class AliRsnCutCentrality : public AliRsnCut {

public:

   AliRsnCutCentrality(const char *name = "cut", const char *est = "VOM", Double_t min = 0, Double_t max = 100.0);
   AliRsnCutCentrality(const AliRsnCutCentrality& copy);
   AliRsnCutCentrality& operator=(const AliRsnCutCentrality& copy);
   virtual ~AliRsnCutCentrality() { }
   
   void           SetEstimator(const char *est) {SetTitle(est);}  
   virtual Bool_t IsSelected(TObject *object);

private:

   ClassDef(AliRsnCutCentrality,1)
};

#endif
