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
#include <TArrayI.h>

#include "AliRsnEventBuffer.h"
#include "AliRsnPair.h"
#include "AliRsnAnalysisTaskSEBase.h"

class AliRsnEvent;

class AliRsnAnalysisSE : public AliRsnAnalysisTaskSEBase
{
  public:
    AliRsnAnalysisSE(const char *name = "AliRsnAnalysisSE", Int_t bufferSize = 1000);
    ~AliRsnAnalysisSE();

    virtual void    InitIOVars();
    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t *option);
    virtual void    Terminate(Option_t *);

    void AddPairMgr(AliRsnPairMgr *pairmgr);
    void AddPairMgrFromConfig(TString configfile);

    void SetBufferSize(const Int_t size) {fBufferSize = size;}
    Int_t GetBufferSize() const {return fBufferSize;}

    TArrayI FindGoodMatches(Int_t iRef, Int_t &found);
    void SetMixingNum(Int_t i) {fMixingNum = i;}
    void SetMixingCut(AliRsnCutSet *cut) {fMixingCut = cut;}


  private:

    AliRsnAnalysisSE(const AliRsnAnalysisSE& copy) :
      AliRsnAnalysisTaskSEBase(copy),
      fDoesMixing(kFALSE),fMixingNum(0),fMixingCut(0x0),fPairMgrs(0),
      fOutList(0x0),fBuffer(0x0),fBufferSize(0) {}
    AliRsnAnalysisSE& operator= (const AliRsnAnalysisSE&) {return *this;}

    void  ProcessEventAnalysis(AliRsnEvent *curEvent);
    void  PostEventProcess(const Short_t &index = 0);

    Bool_t             fDoesMixing;     // flag set to kTRUE if the task contains pairs for mixing
    Int_t              fMixingNum;      // number of events to mix
    AliRsnCutSet      *fMixingCut;      // mixing cut
    TObjArray          fPairMgrs;       // collections of pairs
    TList             *fOutList;        // List of output
    AliRsnEventBuffer *fBuffer;         // event buffer
    Int_t              fBufferSize;     // number of events in buffer

    ClassDef(AliRsnAnalysisSE, 1)
};

#endif
