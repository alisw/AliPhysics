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

#include "AliRsnEvent.h"
#include "AliRsnPIDIndex.h"
#include "AliRsnVATProcessInfo.h"

class AliESDEvent;
class AliAODEvent;
class AliMCEvent;

class AliRsnVAnalysisTaskSE : public AliAnalysisTaskSE
{
  public:
    enum {
      kMaxNumberOfOutputs=10
    };

    AliRsnVAnalysisTaskSE(const char *name = "AliRsnVAnalysisTaskSE", Int_t numOfOutputs = 1, Bool_t mcOnly = kFALSE);
    AliRsnVAnalysisTaskSE(const AliRsnVAnalysisTaskSE& copy);
    AliRsnVAnalysisTaskSE& operator= (const AliRsnVAnalysisTaskSE& /*copy*/) { return *this; }
    virtual ~AliRsnVAnalysisTaskSE() {/* Does nothing*/;}

    virtual void    LocalInit();
    virtual Bool_t  Notify();
    virtual void    ConnectInputData(Option_t *opt);
    // Implementation of interface methods
    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t* opt);
    virtual void    Terminate(Option_t* opt);

    // Implement this
    virtual void    RsnUserCreateOutputObjects();
    virtual void    RsnUserExec(Option_t*);
    virtual void    RsnTerminate(Option_t*);

    virtual void    FillInfo();

    // Prior probs
    AliRsnPIDIndex* GetPIDIndex() {return &fRsnPIDIndex;}
    AliRsnEvent*    GetRsnEvent() {return &fRsnEvent;}
    void            SetPriorProbability(AliPID::EParticleType type, Double_t p) {fRsnEvent.SetPriorProbability(type, p);}
    void            DumpPriors() {fRsnEvent.DumpPriors();}
    void            GetPriorProbability(Double_t *out) const {fRsnEvent.GetPriorProbability(out);}

    void SetMCOnly(Bool_t mcOnly = kTRUE) {fMCOnly = mcOnly;}
    void SetLogType(AliLog::EType_t type, TString otherClasses = "");
    void SetPrintInfoNumber(const Long64_t &num = 100) { fTaskInfo.SetPrintInfoNumber(num); }

  protected:

    AliLog::EType_t         fLogType;
    TString                 fLogClassesString;

    AliESDEvent            *fESDEvent;        // ESD event
    AliMCEvent             *fMCEvent;         // MC event
    AliAODEvent            *fAODEventIn;      // AOD event from input
    AliAODEvent            *fAODEventOut;     // AOD event from output from previous taks

    Bool_t                  fMCOnly;          // use only MC information
    AliRsnEvent             fRsnEvent;        // interface to event for RSN package
    AliRsnPIDIndex          fRsnPIDIndex;     // PID method sorter

    Int_t                   fNumberOfOutputs;
    TList                  *fOutList[kMaxNumberOfOutputs+1]; //!
    AliRsnVATProcessInfo    fTaskInfo;

    void SetDebugForOtherClasses();

    ClassDef(AliRsnVAnalysisTaskSE, 1)
};

#endif
