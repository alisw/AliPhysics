//
// Class AliRsnLoopEffPair
//
// Inherits from basic AliRsnAnalysisTaskEff for efficiency,
// and computed efficiencies for single-tracks
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNLOOPEFFPAIR_H
#define ALIRSNLOOPEFFPAIR_H

#include "AliRsnMother.h"
#include "AliRsnDaughter.h"
#include "AliRsnLoopEff.h"

class AliRsnPairDef;

class AliRsnLoopEffPair : public AliRsnLoopEff {

public:

   AliRsnLoopEffPair(const char *name, AliRsnPairDef *def);
   AliRsnLoopEffPair(const AliRsnLoopEffPair& copy);
   AliRsnLoopEffPair& operator=(const AliRsnLoopEffPair& copy);
   virtual ~AliRsnLoopEffPair() {;}

   AliRsnPairDef* GetDef()                   {return fDef;}
   void           SetDef(AliRsnPairDef *def) {fDef = def;}
   virtual Int_t  DoLoop(AliRsnEvent *main, AliRsnDaughterSelector *smain = 0, AliRsnEvent *mix = 0, AliRsnDaughterSelector *smix = 0);

protected:

   AliRsnPairDef *fDef;         //  used pair definition
   AliRsnMother   fMother;      //! check object (mother)
   AliRsnDaughter fDaughter[2]; //! check object (daughter)

   ClassDef(AliRsnLoopEffPair, 1)
};

#endif


