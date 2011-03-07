//
// Class AliRsnAnalysisTask
//
// Virtual Class derivated from AliRsnVAnalysisTask which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNANALYSISTASK_H
#define ALIRSNANALYSISTASK_H

#include "AliRsnVAnalysisTask.h"
#include "AliRsnAnalysisManager.h"

#include "AliRsnCutSet.h"

class AliRsnPIDDefESD;

class AliRsnAnalysisTask : public AliRsnVAnalysisTask {
public:

   AliRsnAnalysisTask(const char *name = "AliRsnAnalysisTask", Bool_t useKine = kFALSE);
   AliRsnAnalysisTask(const AliRsnAnalysisTask& copy);
   AliRsnAnalysisTask& operator=(const AliRsnAnalysisTask& copy);
   virtual ~AliRsnAnalysisTask() {}

   virtual void            RsnUserCreateOutputObjects();
   virtual void            RsnUserExec(Option_t*);
   virtual void            RsnUserExecMix(Option_t*);
   virtual void            RsnTerminate(Option_t*);
   virtual Bool_t          RsnEventProcess();

   AliRsnCutSet*           GetEventCuts()                           {return &fEventCuts;}
   AliRsnAnalysisManager*  GetAnalysisManager()                     {return &fRsnAnalysisManager;}
   void                    SetAnalysisManagerName(const char *name) {fRsnAnalysisManager.SetName(name);}

   Double_t                GetZeroEventPercentWarning() const            {return fZeroEventPercentWarning;}
   void                    SetZeroEventPercentWarning(Double_t val = 50) {fZeroEventPercentWarning = val;}
   void                    UseZeroEventWarning(Bool_t b = kTRUE)         {fUseZeroEventWarning = b;}

private:

   AliRsnAnalysisManager   fRsnAnalysisManager;      // analysis main engine
   AliRsnCutSet            fEventCuts;               // event cuts
   TList                  *fOutList;                 // list of output events

   Double_t                fZeroEventPercentWarning; // Percent Number for Zero Event Warning
   Bool_t                  fUseZeroEventWarning;     // flag if Zero Event Warning is used (default is true)

   ClassDef(AliRsnAnalysisTask, 1)
};

#endif
