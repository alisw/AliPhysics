//
// Class AliRsnAnalysisTrackEffSE
//
// Virtual Class derivated from AliRsnVAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNANALYSISTRACKEFFSE_H
#define ALIRSNANALYSISTRACKEFFSE_H

#include <TArrayD.h>

#include "AliRsnVAnalysisTaskSE.h"

class AliPID;

class AliCFContainer;

class AliRsnPairDef;
class AliRsnPIDIndex;
class AliRsnPIDDefESD;
class AliRsnCutSet;
class AliRsnFunctionAxis;

class AliRsnAnalysisManager;
class AliRsnAnalysisTrackEffSE : public AliRsnVAnalysisTaskSE
{

  public:
    AliRsnAnalysisTrackEffSE(const char *name = "AliRsnAnalysisTaskEffSE");
    AliRsnAnalysisTrackEffSE(const AliRsnAnalysisTrackEffSE& copy);
    virtual ~AliRsnAnalysisTrackEffSE() {;};

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    // settings
    void            SetEventCuts(AliRsnCutSet *const cuts) {fEventCuts = cuts;}
    void            AddStepMC(AliRsnCutSet *cuts) {fStepListMC.AddLast(cuts);}
    void            AddStepESD(AliRsnCutSet *cuts) {fStepListESD.AddLast(cuts);}
    void            AddAxis(AliRsnFunctionAxis *axis) {fAxisList.AddLast(axis);}

  private:

    AliRsnAnalysisTrackEffSE& operator=(const AliRsnAnalysisTrackEffSE& /*copy*/) {return *this;}
    void                 ProcessEventMC();
    void                 ProcessEventESD();
    Bool_t               PassedAllCutsMC();
    void                 FillContainer(const TObjArray *stepList, Int_t firstOutStep);

    AliRsnCutSet         *fEventCuts;                     // event cuts
    TObjArray             fStepListMC;                    // list of cut steps with MC
    TObjArray             fStepListESD;                   // list of cut steps with ESD
    TObjArray             fAxisList;                      // list of axes of efficiency plots
    AliCFContainer       *fContainer[AliPID::kSPECIES+1]; // one container per particle type + 1 global
    TArrayD               fVar;                           // list of variables of the container
    AliRsnDaughter        fDaughter;                      // interface to track

    ClassDef(AliRsnAnalysisTrackEffSE, 1)
};

#endif
