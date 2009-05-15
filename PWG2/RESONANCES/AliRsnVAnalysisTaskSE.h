//
// Class AliRsnVAnalysisTaskSE
//
// Virtual Class derivated from AliAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNVANALYSISTASKSE_H
#define ALIRSNVANALYSISTASKSE_H

#include <TH1.h>

#include "AliLog.h"

#include "AliAnalysisTaskSE.h"


#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
// class AliRsnEvent;

#include "AliRsnVATProcessInfo.h"

class AliRsnVAnalysisTaskSE : public AliAnalysisTaskSE
{
  public:

    AliRsnVAnalysisTaskSE(const char *name = "AliRsnVAnalysisTaskSE");
    AliRsnVAnalysisTaskSE(const AliRsnVAnalysisTaskSE& copy);
    AliRsnVAnalysisTaskSE& operator= (const AliRsnVAnalysisTaskSE& /*copy*/) { return *this; }
    virtual ~AliRsnVAnalysisTaskSE() {/* Does nothing*/}

    virtual void    LocalInit();
    virtual Bool_t  Notify();
    virtual void    ConnectInputData(Option_t *);
    // Implementation of interface methods
    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t*);
    virtual void    Terminate(Option_t*);

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    virtual void    FillInfo();

    void SetLogType(AliLog::EType_t type, TString otherClasses = "");
    void SetPrintInfoNumber(const Long64_t &num = 100) { fTaskInfo.SetPrintInfoNumber(num); }

  protected:

    AliLog::EType_t         fLogType;
    TString                 fLogClassesString;

    AliESDEvent            *fESDEvent;        // ESD event
    AliMCEvent             *fMCEvent;         // MC event
    AliAODEvent            *fAODEventIn;      // AOD event from input
    AliAODEvent            *fAODEventOut;     // AOD event from output from previous taks

    TList                  *fOutList;
    AliRsnVATProcessInfo    fTaskInfo;
    
    void SetDebugForOtherClasses();

    ClassDef(AliRsnVAnalysisTaskSE, 1)
};

#endif
