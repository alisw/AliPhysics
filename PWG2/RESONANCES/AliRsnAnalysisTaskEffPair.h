//
// Class AliRsnAnalysisTaskEffPair
//
// Inherits from basic AliRsnAnalysisTaskEff for efficiency,
// and computed efficiencies for pairs
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNANALYSISEFFPAIR_H
#define ALIRSNANALYSISEFFPAIR_H

#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnAnalysisTaskEff.h"

class AliRsnPairDef;

class AliRsnAnalysisTaskEffPair : public AliRsnAnalysisTaskEff {

public:

   AliRsnAnalysisTaskEffPair(const char *name = "AliRsnAnalysisTaskMonitorEffSE");
   AliRsnAnalysisTaskEffPair(const AliRsnAnalysisTaskEffPair& copy);
   AliRsnAnalysisTaskEffPair& operator=(const AliRsnAnalysisTaskEffPair& copy);
   virtual ~AliRsnAnalysisTaskEffPair() {;};

protected:

   virtual void    ProcessEventESD();
   virtual void    ProcessEventAOD();
   virtual Int_t   NGoodSteps();
   virtual void    FillContainer(Bool_t mcList);

   AliRsnDaughter  fDaughter[2];   //! current tracks
   AliRsnMother    fMother;        //! current mother
   AliRsnPairDef  *fDef;           //! current definition

   ClassDef(AliRsnAnalysisTaskEffPair, 1)
};

#endif
