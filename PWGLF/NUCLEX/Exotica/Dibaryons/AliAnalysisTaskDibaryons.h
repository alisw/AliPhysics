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
    void SetAnalysisType          ( const char *analysisType  ) { fAnalysisType      = analysisType;      }
    void SetCollidingSystem       ( Int_t  collidingSystem    ) { fCollidingSystem   = collidingSystem;   }
    void SetSelectedTriggerClass  ( UInt_t triggerType        ) { fkTriggerClass     = triggerType;       }
    void SetFilterBit             ( UInt_t filterBit          ) { fFilterBit         = filterBit;         }
    void SetPileupCut             ( Bool_t pileupCut          ) { fPileupCut         = pileupCut;         }
    void SetPairCleaning          ( Bool_t pairCleaning       ) { fPairCleaning      = pairCleaning;      }
    void SetEventMixing           ( Bool_t eventMixing        ) { fEventMixing       = eventMixing;       }
    void SetNsigmaProton          ( Double_t nsigProton       ) { fNsigProton        = nsigProton;        }
    void SetNsigmaV0Daughter      ( Double_t nsigV0Daughter   ) { fNsigV0Daughter    = nsigV0Daughter;    }
    void SetNsigmaCascDaughter    ( Double_t nsigCascDaughter ) { fNsigCascDaughter  = nsigCascDaughter;  }

    void PairCleaner();
    Double_t relKcalc(TLorentzVector track1, TLorentzVector track2);

  private:
    typedef std::deque<TClonesArray*> EventPool;
    typedef std::vector<std::vector<EventPool>> MixingPool;

    TString                 fAnalysisType;            //  "ESD" or "AOD" analysis type
    Int_t                   fCollidingSystem;         //  "pp", "pPb", or "PbPb" colliding system
    UInt_t                  fkTriggerClass;           //  Trigger selection: kINT7, KHighMultV0, etc
    AliPIDResponse         *fPIDResponse;             //! PID response object

    UInt_t                  fFilterBit;               //  filter bit for AOD track selection
    Bool_t                  fPileupCut;               //  apply out-of-bunch pile-up cuts for daughters of V0s and Cascades
    Bool_t                  fPairCleaning;            //  perform Pair Cleaning
    Bool_t                  fEventMixing;             //  perform Event Mixing
    Double_t                fNsigProton;              //  proton PID nsigma
    Double_t                fNsigV0Daughter;          //  V0 daughter PID nsigma
    Double_t                fNsigCascDaughter;        //  Cascade daughter PID nsigma

    THashList              *fOutput;                  //! User output
    AliAODTrack           **fTrackArray;              //! global track info
    TClonesArray           *fProtonArray;             //! proton candidates
    TClonesArray           *fLambdaArray;             //! Lambda candidates
    TClonesArray           *fXiArray;                 //! Xi candidates
    TClonesArray           *fOmegaArray;              //! Omega candidates
    Int_t                   fTrackBuffSize;           //  size of the track array
    Int_t                   fV0BuffSize;              //  size of the V0 array
    Int_t                   fCascadeBuffSize;         //  size of the cascade array

    MixingPool              fProtonEMpool;            //  event mixing pool for proton candidates
    MixingPool              fLambdaEMpool;            //  event mixing pool for Lambda candidates
    MixingPool              fXiEMpool;                //  event mixing pool for Xi candidates
    MixingPool              fOmegaEMpool;             //  event mixing pool for Omega candidates

    AliAnalysisTaskDibaryons(const AliAnalysisTaskDibaryons&);            // not implemented
    AliAnalysisTaskDibaryons& operator=(const AliAnalysisTaskDibaryons&); // not implemented

    ClassDef(AliAnalysisTaskDibaryons, 9);
};

#endif

