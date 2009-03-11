//
// Class AliRsnAnalysisTaskBase
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNANALYSISTASKBASE_H
#define ALIRSNANALYSISTASKBASE_H

#include <TChain.h>

#include "AliAnalysisTask.h"
#include "AliRsnReader.h"

class AliAnalysisManager;
class AliESDEvent;
class AliAODEvent;
class AliRsnEvent;
class AliMCEvent;

#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"


class AliRsnAnalysisTaskBase : public AliAnalysisTask
{
  public:
    AliRsnAnalysisTaskBase(const char *name = "AliRsnAnalysisTaskBase");
    AliRsnAnalysisTaskBase(const AliRsnAnalysisTaskBase& copy):AliAnalysisTask(copy),
        fNumOfEvents(0),fUseAutoHandler(kTRUE),
        fReader(),fPID(),fAnalysisMgr(0x0) {}
    AliRsnAnalysisTaskBase& operator= (const AliRsnAnalysisTaskBase&) {return *this;}
    virtual ~AliRsnAnalysisTaskBase() {/* Does nothing*/}

    enum EInputType
    {
      kAOD = 0,
      kESD,
      kESDMC,
      kESDTPC,
      kESDMCTPC,
      kMC,
      kRSN,
      kLastIndex
    };

    virtual void   InitIOVars();
    virtual void   LocalInit();
    virtual Bool_t Notify();
    virtual void   ConnectInputData(Option_t *);
    virtual void   CreateOutputObjects() {;}
    virtual void   Exec(Option_t *) {;}
    virtual void   Terminate(Option_t *) {;}

    void SetInputType(EInputType type,AliAnalysisManager* am,Bool_t autohandler=kFALSE, Short_t inputIndex=0);
    EInputType GetInputType(Short_t inputIndex=0) { return fInputType[inputIndex]; }

    TChain* GetChain(const Int_t& index = 0) const { return fChain[index]; }

    AliRsnEvent *GetRSNEvent(Int_t index=0) { return fRSN[index]; }

    void SetAnalysisMgr(AliAnalysisManager* theValue) { fAnalysisMgr = theValue; }
    AliAnalysisManager* GetAnalysisMgr() const { return fAnalysisMgr; }

    AliESDInputHandler* GetESDHandler(const Int_t& theValue=0) const { return fRsnESDEH[theValue]; }
    AliMCEventHandler* GetMCHandler(const Int_t& theValue=0) const { return fRsnMCEH[theValue]; }
    AliAODInputHandler* GetAODHandler(const Int_t& theValue=0) const { return fRsnAODEH[theValue]; }

    AliRsnReader *GetReader() { return &fReader; }
    AliRsnPID *GetPID() { return &fPID;}

  protected:

    Long64_t      fNumOfEvents;       // number of events

    TChain       *fChain[2];          // input chain
    EInputType    fInputType[2];      // input type
    Bool_t        fUseAutoHandler;    // flag if should create handler

    AliRsnEvent   *fRSN[2];           // RSN (internal format) event
    AliESDEvent   *fRsnESD[2];        // ESD event
    AliMCEvent    *fRsnMC[2];         // ESD event
    AliAODEvent   *fRsnAOD[2];        // AOD event

    AliESDInputHandler   *fRsnESDEH[2];  // ESD event handler
    AliMCEventHandler    *fRsnMCEH[2];   // ESD event handler
    AliAODInputHandler   *fRsnAODEH[2];  // AOD event handler

    AliRsnReader  fReader;               // Reader
    AliRsnPID     fPID;                  // PID

    AliAnalysisManager *fAnalysisMgr; // pointer to current AnalysisMgr

    virtual void  UseAutoHandler(const Bool_t& theValue) {fUseAutoHandler = theValue;}

    virtual void  ConnectInputDataByInputType(EInputType type ,Short_t inputIndex=0);
    virtual void  ConnectRSN(Short_t inputIndex);
    virtual void  ConnectESD(Short_t inputIndex);
    virtual void  ConnectESDMC(Short_t inputIndex);
    virtual void  ConnectAOD(Short_t inputIndex);

    virtual AliRsnEvent*  GetRsnEventFromInputType(const Short_t &index=0);
    virtual AliRsnEvent*  GetRsnFromAOD(const Short_t &index=0);
    virtual AliRsnEvent*  GetRsnFromESD(const Short_t &index=0);
    virtual AliRsnEvent*  GetRsnFromESDMC(const Short_t &index=0);
    virtual AliRsnEvent*  GetRsnFromRSN(const Short_t &index=0);
    virtual AliRsnEvent*  GetRsnFromMC(const Short_t &index=0);


    ClassDef(AliRsnAnalysisTaskBase, 1)
};

#endif
