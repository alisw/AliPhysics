#ifndef ALIANALYSISTASKDIBARYONS_H
#define ALIANALYSISTASKDIBARYONS_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskDibaryons : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskDibaryons();
    AliAnalysisTaskDibaryons(const char* name);
    virtual ~AliAnalysisTaskDibaryons();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    AliEventCuts fAliEventCuts;
    void SetAnalysisType          (const char *analysisType                    ) { fAnalysisType    = analysisType;    }
    void SetCollidingSystem       (Int_t  collidingSystem                      ) { fCollidingSystem = collidingSystem; }
    void SetSelectedTriggerClass  (AliVEvent::EOfflineTriggerTypes triggerType ) { fkTriggerClass   = triggerType;     }
    void SetPileupCut             (Bool_t pileupCut                            ) { fPileupCut       = pileupCut;       } 

  private:
    TString                 fAnalysisType;            // "ESD" or "AOD" analysis type
    Int_t                   fCollidingSystem;         // "pp" or "pPb" colliding system
    AliVEvent::EOfflineTriggerTypes fkTriggerClass;   // Trigger selection: kINT7, KHighMultV0, etc
    AliPIDResponse         *fPIDResponse;             //! PID response object

    Bool_t                  fPileupCut;               // apply out-of-bunch pile-up cuts for daughters of V0s and Cascades

    THashList              *fOutput;                  //! User output

    AliAnalysisTaskDibaryons(const AliAnalysisTaskDibaryons&);            // not implemented
    AliAnalysisTaskDibaryons& operator=(const AliAnalysisTaskDibaryons&); // not implemented

    ClassDef(AliAnalysisTaskDibaryons, 2);
};

#endif

