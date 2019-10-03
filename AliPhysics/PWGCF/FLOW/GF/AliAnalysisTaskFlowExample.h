#ifndef AliAnalysisTaskFlowExample_H
#define AliAnalysisTaskFlowExample_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"

class AliAnalysisTaskFlowExample : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskFlowExample();
                                AliAnalysisTaskFlowExample(const char *name);
        virtual                 ~AliAnalysisTaskFlowExample();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        //event and track selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
        void                    SetCentralityEst(TString est){ fCentEstimator = est; }
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        AliEventCuts            fEventCuts;
        void                    SetPtRange(Double_t min, Double_t max) {fPtMin = min; fPtMax = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        //output histograms
        TH2F*                   fHistPhiEta;    //!

        AliAnalysisTaskFlowExample(const AliAnalysisTaskFlowExample&); // not implemented
        AliAnalysisTaskFlowExample& operator=(const AliAnalysisTaskFlowExample&); // not implemented

        //event and track selection
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fEventRejectAddPileUp;
        UInt_t                  fFilterBit;
        Double_t                fPtMin;
        Double_t                fPtMax;
        Double_t                fAbsEtaMax;
        TString                 fCentEstimator;
        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;

        ClassDef(AliAnalysisTaskFlowExample, 1);
};

#endif
