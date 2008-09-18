//
// Class AliRsnAnalysisSE
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNANALYSISAT_H
#define ALIRSNANALYSISAT_H

#include <TH1.h>

#include "AliRsnEventBuffer.h"
#include "AliRsnPair.h"
#include "AliRsnAnalysisTaskSEBase.h"

class AliRsnEvent;

class AliRsnAnalysisSE : public AliRsnAnalysisTaskSEBase
{
  public:
    AliRsnAnalysisSE(const char *name = "AliRsnAnalysisSE");
    AliRsnAnalysisSE(const AliRsnAnalysisSE& copy): AliRsnAnalysisTaskSEBase(copy),
       fPairMgrs(0),fOutList(0x0),fRsnEventBuffer(0x0),fNumOfEventsInBuffer(100) {}
    AliRsnAnalysisSE& operator= (const AliRsnAnalysisSE&) {return *this;}
    ~AliRsnAnalysisSE();

    virtual void    InitIOVars();
//     virtual void    LocalInit();
//     virtual Bool_t  Notify();
    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t *option);
    virtual void    Terminate(Option_t *);

    void AddPairMgr(AliRsnPairMgr*pairmgr);

    void SetNumOfEventsInBuffer(const Int_t& theValue) { fNumOfEventsInBuffer = theValue; }
    Int_t GetNumOfEventsInBuffer() const { return fNumOfEventsInBuffer; }


  private:

    TObjArray         fPairMgrs;

    TList             *fOutList;              // List of output
    AliRsnEventBuffer *fRsnEventBuffer;       // event buffer
    Int_t             fNumOfEventsInBuffer;  // number of events in buffer

    void            ProcessEventAnalysis(AliRsnEvent *curEvent);
    void            PostEventProcess(const Short_t &index=0);

    ClassDef(AliRsnAnalysisSE, 1)
};

#endif
