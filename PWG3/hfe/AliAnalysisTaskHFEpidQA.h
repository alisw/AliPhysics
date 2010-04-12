#ifndef ALIANALYSISTASKHFEPIDQA_H
#define ALIANALYSISTASKHFEPIDQA_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TH1;
class TList;

class AliHFEpidQA;

class AliAnalysisTaskHFEpidQA : public AliAnalysisTaskSE{
  enum{
    kV0pidQA = BIT(22),
    kRecalculateTRDpid = BIT(23)
  };
  public:
    AliAnalysisTaskHFEpidQA();
    AliAnalysisTaskHFEpidQA(const Char_t *name);
    ~AliAnalysisTaskHFEpidQA();

    void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    Bool_t HasV0pidQA() const { return TestBit(kV0pidQA); };
    Bool_t HasRecalculateTRDpid() const { return TestBit(kRecalculateTRDpid); };
    void SetV0pidQA(Bool_t v0pidQA = kTRUE) { SetBit(kV0pidQA, v0pidQA); };
    void SetRecalculateTRDpid(Bool_t recal = kTRUE) { SetBit(kRecalculateTRDpid, recal); };

  private:
    AliAnalysisTaskHFEpidQA(const AliAnalysisTaskHFEpidQA &ref);
    AliAnalysisTaskHFEpidQA &operator=(const AliAnalysisTaskHFEpidQA &ref);
    AliHFEpidQA *fPIDqa;    //! The heart of the analysis  
    TList *fOutput;         //! Container for output histos
    TH1 *fEvents;           //! Number of Events

    ClassDef(AliAnalysisTaskHFEpidQA, 1)
};

#endif

