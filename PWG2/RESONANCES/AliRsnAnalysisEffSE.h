//
// Class AliRsnAnalysisEffSE
//
// Virtual Class derivated from AliRsnVAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNANALYSISEFFSE_H
#define ALIRSNANALYSISEFFSE_H

#include <TArrayD.h>

#include "AliRsnVAnalysisTaskSE.h"
#include "AliRsnEvent.h"
#include "AliRsnMother.h"

class AliPID;

class AliCFContainer;

class AliRsnPairDef;
class AliRsnPIDIndex;
class AliRsnPIDDefESD;
class AliRsnCutSet;
class AliRsnCutManager;
class AliRsnValue;

class AliRsnAnalysisManager;
class AliRsnAnalysisEffSE : public AliRsnVAnalysisTaskSE
{

  public:
    AliRsnAnalysisEffSE(const char *name = "AliRsnAnalysisTaskEffSE");
    AliRsnAnalysisEffSE(const AliRsnAnalysisEffSE& copy);
    virtual ~AliRsnAnalysisEffSE() {;};

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    // settings
    void            SetEventCuts(AliRsnCutSet *const cuts) {fEventCuts = cuts;}
    void            AddPairDef(AliRsnPairDef *pairDef);
    void            AddStepMC(AliRsnCutManager *mgr)  {fStepListMC.AddLast(mgr);}
    void            AddStepESD(AliRsnCutManager *mgr) {fStepListESD.AddLast(mgr);}
    void            AddAxis(AliRsnValue *axis) {fAxisList.AddLast(axis);}

  private:

    AliRsnAnalysisEffSE& operator=(const AliRsnAnalysisEffSE& /*copy*/) {return *this;}
    void                 ProcessEvent(AliRsnPairDef *pairDef);
    void                 ProcessEventMC(AliRsnPairDef *pairDef);
    void                 ProcessEventESD(AliRsnPairDef *pairDef);
    void                 FillContainer(AliCFContainer *cont, const TObjArray *stepList, AliRsnPairDef *pd, Int_t firstOutStep);
    Int_t                FindESDtrack(Int_t label, AliESDEvent *esd, Bool_t rejectFakes);

    Bool_t                fUseITSSA;                // switch to use ITS standalone tracks
    Bool_t                fUseGlobal;               // switch to use global tracks
    AliRsnCutSet         *fEventCuts;               // event cuts
    TObjArray             fStepListMC;              // list of cut managers for all steps with MC
    TObjArray             fStepListESD;             // list of cut managers for all steps with ESD
    TObjArray             fAxisList;                // list of axes of efficiency plots
    TObjArray             fPairDefList;             // decay channels
    TList                *fContainerList;           // list of CF containers
    TList                *fOutList;                 // global output list
    TArrayD               fVar;                     // list of variables of the container
    AliRsnMother          fPair;                    // interface to pair
    AliRsnDaughter        fDaughter[2];             // interface to tracks

    ClassDef(AliRsnAnalysisEffSE, 1)
};

#endif
