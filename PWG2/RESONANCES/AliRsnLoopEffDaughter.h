//
// Class AliRsnLoopEffDaughter
//
// Inherits from basic AliRsnAnalysisTaskEff for efficiency,
// and computed efficiencies for single-tracks
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNLOOPEFFDAUGHTER_H
#define ALIRSNLOOPEFFDAUGHTER_H

#include "AliRsnLoopEff.h"

class AliRsnDaughterDef;

class AliRsnLoopEffDaughter : public AliRsnLoopEff {

public:

   AliRsnLoopEffDaughter(const char *name, AliRsnDaughterDef *def);
   AliRsnLoopEffDaughter(const AliRsnLoopEffDaughter& copy);
   AliRsnLoopEffDaughter& operator=(const AliRsnLoopEffDaughter& copy);
   virtual ~AliRsnLoopEffDaughter() {;};

   AliRsnDaughterDef* GetDef()                       {return fDef;}
   void               SetDef(AliRsnDaughterDef *def) {fDef = def;}

   virtual Bool_t     OkStepMC(TObject *checked, Int_t step);
   virtual Bool_t     OkStepRec(TObject *checked, Int_t step);

protected:

   virtual Int_t ProcessEventESD(AliRsnEvent *rsn);
   virtual Int_t ProcessEventAOD(AliRsnEvent *rsn);

   AliRsnDaughterDef *fDef;  // used daughter definition

   ClassDef(AliRsnLoopEffDaughter, 1)
};

#endif

