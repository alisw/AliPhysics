#ifndef AliAnalysisTaskCorrForFlow_H
#define AliAnalysisTaskCorrForFlow_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliLog.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TList.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TString.h"
#include <vector>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliEventPoolManager.h"
#include "AliVEvent.h"
#include "AliTHn.h"

class AliAnalysisTaskCorrForFlow : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskCorrForFlow();
                                AliAnalysisTaskCorrForFlow(const char *name);
        virtual                 ~AliAnalysisTaskCorrForFlow();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        //event and track selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetPtRangeTrig(Double_t min, Double_t max) {fPtMinTrig = min; fPtMaxTrig = max; }
        void                    SetPtRangeAss(Double_t min, Double_t max) {fPtMinAss = min; fPtMaxAss = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetCentrality(TString cent, Double_t min = 0.0, Double_t max = 20.0) { fCentEstimator = cent; fCentMin = min; fCentMax = max; }
        void                    SetPtBins(std::vector<Double_t> bins) { fPtBinsTrigCharged = bins; }
        void                    SetPtBinsAss(std::vector<Double_t> bins) { fPtBinsAss = bins; }

    private:

        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;

        Double_t                RangePhi(Double_t DPhi);
        Double_t                GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius);
        void                    FillCorrelations();
        void                    FillCorrelationsMixed();


        AliAnalysisTaskCorrForFlow(const AliAnalysisTaskCorrForFlow&); // not implemented
        AliAnalysisTaskCorrForFlow& operator=(const AliAnalysisTaskCorrForFlow&); // not implemented

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputListCharged;    //! output list
        TObjArray*              fTracksTrigCharged; //!
        TObjArray*              fTracksAss; //!
        AliEventPoolManager*    fPoolMgr;  //!  event pool manager for Event Mixing
        //output histograms
        TH1D*                   fhEventCounter; //!
        TH2D*                   fHistPhiEta; //!
        TH2D*                   fhTrigTracks; //!
        AliTHn*                 fhChargedSE; //!
        AliTHn*                 fhChargedME; //!


        //event and track selection
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fIsHMpp;
        UInt_t                  fFilterBit;
        Int_t                   fbSign;
        Double_t                fPtMinTrig;
        Double_t                fPtMaxTrig;
        Double_t                fPtMinAss;
        Double_t                fPtMaxAss;
        std::vector<Double_t>   fPtBinsTrigCharged;
        std::vector<Double_t>   fPtBinsAss;
        Double_t                fCentMin;
        Double_t                fCentMax;
        Double_t                fCentrality;
        Double_t                fAbsEtaMax;
        Double_t                fPVz;
        TString                 fCentEstimator;
        AliEventCuts            fEventCuts;

        // mixing
        Int_t                   fPoolMaxNEvents;   // maximum number of events in the pool
        Int_t                   fPoolMinNTracks;   // minimum number of tracks in the pool
        Int_t                   fMinEventsToMix;   // minimum number of events for mixing
        Int_t                   fNzVtxBins; // number of PV z bins
        Int_t                   fNCentBins; // number of centrality bins
        std::vector<Double_t>   fzVtxBins;
        std::vector<Double_t>   fCentBins;
        Double_t                fMergingCut; // [0.02] cut for track spliting/merging

        ClassDef(AliAnalysisTaskCorrForFlow, 1);
};

#endif
