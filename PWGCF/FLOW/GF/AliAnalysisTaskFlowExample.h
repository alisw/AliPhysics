#ifndef AliAnalysisTaskFlowExample_H
#define AliAnalysisTaskFlowExample_H

#include <vector>
#include <complex>

#include "TChain.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TFile.h"
#include "TProfile.h"

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"

class AliAnalysisTaskFlowExample : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskFlowExample();
                                AliAnalysisTaskFlowExample(const char *name, Bool_t bUseWeights);
        virtual                 ~AliAnalysisTaskFlowExample();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        //event selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
        void                    SetCentralityEst(TString est){ fCentEstimator = est; }
        void                    SetCentralityBin(Int_t nOfBins, Double_t cenMin, Double_t cenMax){ fCentBinNum = nOfBins; fCentMin = cenMin; fCentMax = cenMax; }
        //track selection
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetPtRange(Double_t min, Double_t max) {fPtMin = min; fPtMax = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetEtaGap(Double_t eta) {fEtaGap = eta; }
        void                    SetHarmonic(Int_t harm) {fHarmonic = harm; }
        void                    SetIsHMPP(Bool_t ispp = kTRUE) {fIsPP = ispp; }


    private:
        AliAnalysisTaskFlowExample(const AliAnalysisTaskFlowExample&); // not implemented
        AliAnalysisTaskFlowExample& operator=(const AliAnalysisTaskFlowExample&); // not implemented

        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;

        //flow ana
        static const Int_t      maxHarm = 7;
        static const Int_t      maxPow = 7;
        void                    ResetQvector(std::complex<double> (&array)[maxHarm][maxPow]);
        void                    FillQvector(std::complex<double> (&array)[maxHarm][maxPow], const Double_t dPhi, const Double_t dWeight);
        void                    CalculateCorrelations();
        Bool_t                  AreWeightsLoaded();
        Double_t                GetWeight(const Double_t dPhi);

        std::complex<double>     Q(Int_t n, Int_t p) const;
        std::complex<double>     QGapPos(Int_t n, Int_t p) const;
        std::complex<double>     QGapNeg(Int_t n, Int_t p) const;
        std::complex<double>     Two(Int_t n1, Int_t n2) const;
        std::complex<double>     TwoGap(Int_t n1, Int_t n2) const;
        std::complex<double>     Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4) const;
        std::complex<double>     Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6) const;


        AliAODEvent*            fAOD;           //! input event
        //event and track selection
        AliVEvent::EOfflineTriggerTypes    fTrigger; // [AliVEvent::kINT7]
        Bool_t                  fIsPP; // [kFALSE]
        Bool_t                  fEventRejectAddPileUp; // [kFALSE]
        Bool_t                  fUseWeights; // [kFALSE]
        UInt_t                  fFilterBit; // [96]
        Int_t                   fCentBinNum; // [10]
        Double_t                fCentMin; // [0.0]
        Double_t                fCentMax; // [100.0]
        Double_t                fCentrality; // [1]
        Double_t                fPtMin; // [0.2]
        Double_t                fPtMax; // [5.0]
        Double_t                fAbsEtaMax; // [1.0]
        Double_t                fPVtxCutZ; // [10.0] PV z cut (in cm)
        Double_t                fPVz; // [99.9] PV z-coordinate
        Double_t                fEtaGap; // [-1.0]
        TString                 fCentEstimator; // [kV0M]
        AliEventCuts            fEventCuts;

        // flow ana
        std::complex<double>    Qvector[maxHarm][maxPow];
        std::complex<double>    QvectorPos[maxHarm][maxPow];
        std::complex<double>    QvectorNeg[maxHarm][maxPow];
        Int_t                   fHarmonic; // [2]

        //outputs
        TList*                  fOutputList;    //! output list
        TList*                  fInputWeights;    //! list w weights
        //output events QA
        TH1D*                   fhEventCounter; //!
        TH1D*                   fhWeight; //!
        TH2F*                   fHistPhiEta;    //!
        TProfile*               fProf; //!
        TProfile*               fProfGap; //!
        TProfile*               fProf4PC; //!
        TProfile*               fProf6PC; //!

        ClassDef(AliAnalysisTaskFlowExample, 2);
};

#endif
