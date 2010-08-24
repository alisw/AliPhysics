//
// Class AliRsnAnalysisSE
//
// Virtual Class derivated from AliRsnVAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNANALYSISSE_H
#define ALIRSNANALYSISSE_H

#include "AliPID.h"
#include "AliRsnVAnalysisTaskSE.h"
#include "AliRsnAnalysisManager.h"
#include "AliRsnEvent.h"
#include "AliRsnCutSet.h"

class AliRsnPIDDefESD;

class AliRsnAnalysisSE : public AliRsnVAnalysisTaskSE
{

  public:
    AliRsnAnalysisSE(const char *name = "AliRsnAnalysisSE", Bool_t useKine = kFALSE);
    AliRsnAnalysisSE(const AliRsnAnalysisSE& copy);
//     virtual ~AliRsnAnalysisSE();

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    AliRsnAnalysisManager *GetAnalysisManager() {return &fRsnAnalysisManager;}
    void                   SetAnalysisManagerName(const char *name) {fRsnAnalysisManager.SetName(name);}

    AliRsnCutSet* GetEventCuts() {return &fEventCuts;}
//     void          SetEventCuts(AliRsnCutSet *const cuts) {fEventCuts = cuts;}

    Double_t GetZeroEventPercentWarning() const {return fZeroEventPercentWarning;}
    void     SetZeroEventPercentWarning(Double_t val = 50) {fZeroEventPercentWarning = val;}
    void     UseZeroEventWarning(Bool_t b = true) {fUseZeroEventWarning = b;}

  private:

    AliRsnAnalysisSE& operator=(const AliRsnAnalysisSE& /*copy*/) {return *this;}

    AliRsnAnalysisManager fRsnAnalysisManager;  // analysis main engine
    AliRsnCutSet          fEventCuts;           // event cuts
    TList                *fOutList;             // list of output events

    Double_t              fZeroEventPercentWarning; // Percent Number for Zero Event Warning
    Bool_t                fUseZeroEventWarning;     // flag if Zero Event Warning is used (default is true)

    ClassDef(AliRsnAnalysisSE, 1)
};

#endif
