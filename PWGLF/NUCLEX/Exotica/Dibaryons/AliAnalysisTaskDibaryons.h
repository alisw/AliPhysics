#ifndef ALIANALYSISTASKDIBARYONS_H
#define ALIANALYSISTASKDIBARYONS_H

#include <vector>
#include <deque>
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
    void SetFilterBit             (UInt_t filterBit                            ) { fFilterBit       = filterBit;       }
    void SetPileupCut             (Bool_t pileupCut                            ) { fPileupCut       = pileupCut;       } 
    void SetEventMixing           (Bool_t eventMixing                          ) { fEventMixing     = eventMixing;     }

    void PairCleaner();
    Double_t relKcalc(TLorentzVector track1, TLorentzVector track2);

  private:
    typedef std::deque<TClonesArray*> EventPool;
    typedef std::vector<std::vector<EventPool>> MixingPool;

    TString                 fAnalysisType;            // "ESD" or "AOD" analysis type
    Int_t                   fCollidingSystem;         // "pp" or "pPb" colliding system
    AliVEvent::EOfflineTriggerTypes fkTriggerClass;   // Trigger selection: kINT7, KHighMultV0, etc
    AliPIDResponse         *fPIDResponse;             //! PID response object

    UInt_t                  fFilterBit;               // filter bit for AOD track selection
    Bool_t                  fPileupCut;               // apply out-of-bunch pile-up cuts for daughters of V0s and Cascades
    Bool_t                  fEventMixing;             //

    THashList              *fOutput;                  //! User output
    AliAODTrack           **fTrackArray;              //! global track info
    TClonesArray           *fProtonArray;             //! proton candidates
    TClonesArray           *fLambdaArray;             //! Lambda candidates
    TClonesArray           *fXiArray;                 //! Xi candidates
    TClonesArray           *fOmegaArray;              //! Omega candidates
    const Int_t             fTrackBuffSize;           // size of the track array
    const Int_t             fV0BuffSize;              // size of the V0 array
    const Int_t             fCascadeBuffSize;         // size of the cascade array

    MixingPool              fProtonEMpool;            //
    MixingPool              fLambdaEMpool;            //
    MixingPool              fXiEMpool;                //
    MixingPool              fOmegaEMpool;             //

    AliAnalysisTaskDibaryons(const AliAnalysisTaskDibaryons&);            // not implemented
    AliAnalysisTaskDibaryons& operator=(const AliAnalysisTaskDibaryons&); // not implemented

    ClassDef(AliAnalysisTaskDibaryons, 7);
};

#endif

