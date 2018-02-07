#ifndef AliAodSkimTask_H
#define AliAodSkimTask_H

#include <AliAnalysisTaskSE.h>
class AliAODMCHeader;
class TH1F;

class AliAodSkimTask: public AliAnalysisTaskSE  
{
  public:
    AliAodSkimTask();
    AliAodSkimTask(const char *name);
    virtual              ~AliAodSkimTask();
    void                  SetClusMinE(Double_t v) {fClusMinE=v;}
    void                  SetCutMC(Bool_t b)      {fCutMC=b;   }
  protected:
    void                  UserCreateOutputObjects();
    void                  UserExec(Option_t* option);
    void                  Terminate(Option_t* option);
    Double_t              fClusMinE;      //  minimum cluster energy to accept event
    Bool_t                fCutMC;         //  if true cut MC particles with |Y|>1.2
    UInt_t                fTrials;        //! events seen since last acceptance 
    AliAODEvent          *fAOD;           //! input event
    AliAODMCHeader       *fAODMcHeader;   //! MC header
    TList                *fOutputList;    //! output list
    TH1F                 *fHevs;          //! events processed/accepted
    TH1F                 *fHclus;         //! cluster distribution
    AliAodSkimTask(const AliAodSkimTask&);             // not implemented
    AliAodSkimTask& operator=(const AliAodSkimTask&);  // not implemented
  ClassDef(AliAodSkimTask, 1);
};
#endif
