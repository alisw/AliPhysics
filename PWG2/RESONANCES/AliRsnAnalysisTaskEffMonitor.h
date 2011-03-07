//
// Class AliRsnAnalysisTaskEffMonitor
//
// Inherits from basic AliRsnAnalysisTaskEff for efficiency,
// and computed efficiencies for single-tracks
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNANALYSISEFFMONITOR_H
#define ALIRSNANALYSISEFFMONITOR_H

#include "AliRsnDaughter.h"
#include "AliRsnAnalysisTaskEff.h"

class AliRsnDaughterDef;

class AliRsnAnalysisTaskEffMonitor : public AliRsnAnalysisTaskEff {

public:

   AliRsnAnalysisTaskEffMonitor(const char *name = "AliRsnAnalysisTaskMonitorEffSE");
   AliRsnAnalysisTaskEffMonitor(const AliRsnAnalysisTaskEffMonitor& copy);
   AliRsnAnalysisTaskEffMonitor& operator=(const AliRsnAnalysisTaskEffMonitor& copy);
   virtual ~AliRsnAnalysisTaskEffMonitor() {;};

protected:

   virtual void    ProcessEventESD();
   virtual void    ProcessEventAOD();
   virtual Int_t   NGoodSteps();
   virtual void    FillContainer(Bool_t mcList, TObject *def);

   AliRsnDaughter     fDaughter;   //! current track

   ClassDef(AliRsnAnalysisTaskEffMonitor, 1)
};

#endif
