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

#include "AliAnalysisTaskSE.h"

#include "AliRsnEvent.h"
#include "AliRsnVATProcessInfo.h"

class AliESDEvent;
class AliAODEvent;
class AliMCEvent;

class AliRsnVAnalysisTaskSE : public AliAnalysisTaskSE
{
  public:
  
    AliRsnVAnalysisTaskSE(const char *name = "AliRsnVAnalysisTaskSE", Bool_t mcOnly = kFALSE);
    AliRsnVAnalysisTaskSE(const AliRsnVAnalysisTaskSE& copy);
    AliRsnVAnalysisTaskSE& operator= (const AliRsnVAnalysisTaskSE& /*copy*/) { return *this; }
    virtual ~AliRsnVAnalysisTaskSE() {/* Does nothing*/;}

    // basic interface methods
    virtual void    LocalInit();
    virtual Bool_t  UserNotify();
    virtual void    ConnectInputData(Option_t *opt);
    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t* opt);
    virtual void    Terminate(Option_t* opt);

    // customized methods (to be implemented in derived classes)
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    // event pre-processing functions
    virtual Bool_t  EventProcess();

    // getters
    AliRsnEvent*           GetRsnEvent() {return &fRsnEvent;}
    AliRsnVATProcessInfo*  GetInfo()     {return &fTaskInfo;}

    // setters
    void SetMCOnly(Bool_t mcOnly = kTRUE)                           {fMCOnly = mcOnly;}
    void SetLogType(AliLog::EType_t type, const char *classes = "") {fLogType = type; fLogClassesString = classes;}
    void SetPrintInfoNumber(const Long64_t &num = 100)              {fTaskInfo.SetPrintInfoNumber(num);}

  protected:

    AliLog::EType_t         fLogType;          //  log type
    TString                 fLogClassesString; //  all classes string divided with ":"

    AliESDEvent            *fESDEvent;         //  ESD event
    AliMCEvent             *fMCEvent;          //  MC event
    AliAODEvent            *fAODEventIn;       //  AOD event from input
    AliAODEvent            *fAODEventOut;      //  AOD event from output from previous taks

    Bool_t                  fMCOnly;           //  use only MC information
    AliRsnEvent             fRsnEvent;         //  interface to event for RSN package

    TList                  *fInfoList;         //! output list for informations
    AliRsnVATProcessInfo    fTaskInfo;         //  task info

    void                    SetDebugForAllClasses();

    ClassDef(AliRsnVAnalysisTaskSE, 1)
};

#endif
