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
#include "AliRsnPIDIndex.h"
#include "AliRsnEvent.h"

class AliRsnPIDDefESD;
class AliRsnCutSet;

class AliRsnAnalysisSE : public AliRsnVAnalysisTaskSE
{

  public:
    AliRsnAnalysisSE(const char *name = "AliRsnAnalysisSE",Int_t numOfOutputs=1,Bool_t useKine=kFALSE);
    AliRsnAnalysisSE(const AliRsnAnalysisSE& copy);
    virtual ~AliRsnAnalysisSE() {;};

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    AliRsnAnalysisManager *GetAnalysisManager(Int_t index = 0, TString name = "");
    void                   SetAnalysisManagerName(const char *name, Int_t index = 0)
                            {fRsnAnalysisManager[index].SetName(name);}

    AliRsnCutSet* GetEventCuts() const {return fEventCuts;}
    void          SetEventCuts(AliRsnCutSet *const cuts) {fEventCuts = cuts;}

    Double_t GetZeroEventPercentWarning() const { return fZeroEventPercentWarning;}
    void     SetZeroEventPercentWarning(Double_t val = 50) { fZeroEventPercentWarning = val;}
    void     UseZeroEventWarning(Bool_t b = true) { fUseZeroEventWarning = b;}

  private:

    AliRsnAnalysisSE& operator=(const AliRsnAnalysisSE& /*copy*/) {return *this;}

    AliRsnAnalysisManager fRsnAnalysisManager[10];  // analysis main engine
    AliRsnCutSet         *fEventCuts;               // event cuts

    Double_t              fZeroEventPercentWarning; //! Percent Number for Zero Event Warning
    Bool_t                fUseZeroEventWarning;     //! flag if Zero Event Warning is used (default is true)

    ClassDef(AliRsnAnalysisSE, 1)
};

#endif
