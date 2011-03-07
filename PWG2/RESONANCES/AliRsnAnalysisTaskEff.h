//
// Class AliRsnAnalysisTaskEff
//
// Base class for efficiency computation tasks
// which should be inherited by different efficiency computators
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNANALYSISTASKEFF_H
#define ALIRSNANALYSISTASKEFF_H

#include <TArrayI.h>
#include <TArrayD.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliCFContainer.h"

#include "AliRsnValue.h"
#include "AliRsnCutSet.h"
#include "AliRsnVAnalysisTask.h"

class TList;
class AliVEvent;

class AliRsnAnalysisTaskEff : public AliRsnVAnalysisTask {

public:

   AliRsnAnalysisTaskEff(const char *name = "AliRsnAnalysisTasEff");
   AliRsnAnalysisTaskEff(const AliRsnAnalysisTaskEff& copy);
   AliRsnAnalysisTaskEff& operator=(const AliRsnAnalysisTaskEff& copy);
   virtual ~AliRsnAnalysisTaskEff() {;};
   
   // work-flow
   AliRsnCutSet*   GetEventCuts() {return &fEventCuts;}
   void            AddDef(TObject *def);
   void            AddAxis(AliRsnValue *axis);
   void            AddStepMC(TObject *set);
   void            AddStepRec(TObject *set);

   // inherited
   virtual void    RsnUserCreateOutputObjects();
   virtual void    RsnUserExec(Option_t*);
   virtual void    RsnTerminate(Option_t*);
   virtual Bool_t  RsnEventProcess();

protected:

   TArrayI         FindTracks(Int_t label, AliVEvent *esd);
   virtual void    ProcessEventESD();
   virtual void    ProcessEventAOD();
   virtual Int_t   NGoodSteps();
   virtual void    FillContainer(Bool_t mcList, TObject *def);
   
   TObjArray       fDefs;        //  list of definitions
   TObjArray       fStepsMC;     //  list of cuts for all steps with MC tracks
   TObjArray       fStepsRec;    //  list of cuts for all steps with reconstructed tracks
   TClonesArray    fAxes;        //  list of axes of efficiency plots
   
   TList          *fOutList;     //  global output list
   AliRsnCutSet    fEventCuts;   //  event cuts
   
   TArrayD         fVar;         //! list of variables of the container (temporary)

   ClassDef(AliRsnAnalysisTaskEff, 1)
};

#endif
