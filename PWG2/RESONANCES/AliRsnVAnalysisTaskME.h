//
// Class AliRsnVAnalysisTaskME
//
// Virtual Class derivated from AliAnalysisTaskME which will be base class
// for all RSN ME tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef AliRsnVATME_H
#define AliRsnVATME_H

#include <TH1.h>

#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

#include "AliRsnVATProcessInfo.h"

#include "AliAnalysisTaskME.h"

class AliRsnVAnalysisTaskME : public AliAnalysisTaskME
{
public:
    AliRsnVAnalysisTaskME(const char *name = "AliRsnVAnalysisTaskME");
    AliRsnVAnalysisTaskME(const AliRsnVAnalysisTaskME& copy);
    AliRsnVAnalysisTaskME& operator= (const AliRsnVAnalysisTaskME& /*copy*/) {
        return *this;
    }
    virtual ~AliRsnVAnalysisTaskME() {/* Does nothing*/}

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

    void            SetLogType(AliLog::EType_t type,TString otherClasses="");
    void            SetPrintInfoNumber(const Long64_t &num=100) { fTaskInfo.SetPrintInfoNumber(num); }

protected:

    AliLog::EType_t         fLogType;
    TString                 fLogClassesString;

    AliESDEvent            *fESDEvent;        // AliVEvent event
    AliMCEvent             *fMCEvent;         // ESD event
    AliAODEvent            *fAODEvent;        // AOD event

    TList                  *fOutList;
    AliRsnVATProcessInfo    fTaskInfo;
    
    void            SetDebugForOtherClasses();

    ClassDef(AliRsnVAnalysisTaskME, 1)
};

#endif
