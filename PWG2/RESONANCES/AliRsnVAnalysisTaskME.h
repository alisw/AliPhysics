//
// Class AliRsnVAnalysisTaskME
//
// Virtual Class derivated from AliAnalysisTaskME which will be base class
// for all RSN ME tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNVANALYSISTASKME_H
#define ALIRSNVANALYSISTASKME_H

#include "AliLog.h"
#include "AliRsnVATProcessInfo.h"

#include "AliAnalysisTaskME.h"

class TH1;
class AliESDEvent;
class AliMCEvent;
class AliAODEvent;
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
    virtual void    ConnectInputData(Option_t *opt);
    // Implementation of interface methods
    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t *opt);
    virtual void    Terminate(Option_t *opt);

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t *opt);
    virtual void    RsnTerminate(Option_t *opt);

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
