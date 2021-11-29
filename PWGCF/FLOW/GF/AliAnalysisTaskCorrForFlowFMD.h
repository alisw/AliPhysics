#ifndef ALIANALYSISTASKCORRFORFLOWFMD_H
#define ALIANALYSISTASKCORRFORFLOWFMD_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliLog.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TList.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TString.h"
#include "TRandom.h"
#include <vector>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliEventPoolManager.h"
#include "AliVEvent.h"
#include "AliTHn.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODForwardMult.h"
#include "AliAODVZERO.h"
#include "AliPartSimpleForCorr.h"

class AliAnalysisTaskCorrForFlowFMD : public AliAnalysisTaskSE
{
    public:

                                AliAnalysisTaskCorrForFlowFMD();
                                AliAnalysisTaskCorrForFlowFMD(const char *name, Bool_t bUseEff);
        virtual                 ~AliAnalysisTaskCorrForFlowFMD();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        enum                    AnaType { eTPCFMDA = 0, eTPCFMDC, eFMDAFMDC, eTPCTPC };
        enum                    PartSpecies { eCharged = 0, ePion, eKaon, eProton };

        //event and track selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetPtRangeTrig(Double_t min, Double_t max) {fPtMinTrig = min; fPtMaxTrig = max; }
        void                    SetPtRangeAss(Double_t min, Double_t max) {fPtMinAss = min; fPtMaxAss = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetPVZcut(Double_t cut) {fPVzCut = cut; }
        void                    SetBayesPIDcut(Double_t pion, Double_t kaon, Double_t proton){ fPIDbayesPion = pion; fPIDbayesKaon = kaon; fPIDbayesProton = proton; }
        void                    SetPhiStarCur(Double_t phiStar) {fMergingCut = phiStar; }
        void                    SetCentrality(TString cent, Double_t min = 0.0, Double_t max = 20.0) { fCentEstimator = cent; fCentMin = min; fCentMax = max; }
        void                    SetUseNchRange(Bool_t range, Int_t min, Int_t max) { fUseNch = range; fNchMin = min; fNchMax = max; }
        void                    SetPtBins(std::vector<Double_t> bins) { fPtBinsTrigCharged = bins; }
        void                    SetPtBinsAss(std::vector<Double_t> bins) { fPtBinsAss = bins; }
        void                    SetCentBinsForMixing(Int_t nofBins, std::vector<Double_t> bins) { fNCentBins = nofBins; fCentBins = bins; }
        void                    SetDoPID(Bool_t pid = kTRUE) { fDoPID = pid; }
        void                    SetIsHMpp(Bool_t hm = kTRUE) { fIsHMpp = hm; }
        void                    SetNofSamples(Int_t n) { fNOfSamples = n; }
        void                    SetUseEtaDependentEfficiencies(Bool_t ef = kTRUE) { fEfficiencyEtaDependent = ef; }
        void                    SetUseOppositeSidesOnly(Bool_t sides = kTRUE) { fUseOppositeSidesOnly = sides; }

        // FMD related
        void                    SetAnaType(AnaType type) { fAnalType = type; }
        void                    SetUseFMDcut(Bool_t cut = kTRUE) { fUseFMDcut = cut; }
        void                    SetFMDcutParameters(Double_t par0a, Double_t par1a, Double_t par0c, Double_t par1c) { fFMDcutapar0 = par0a; fFMDcutapar1 = par1a; fFMDcutcpar0 = par0c; fFMDcutcpar1 = par1c; }
        void                    SetFMDacceptanceCuts(Double_t cutAlower, Double_t cutAupper, Double_t cutClower, Double_t cutCupper) { fFMDAacceptanceCutLower = cutAlower; fFMDAacceptanceCutUpper = cutAupper; fFMDCacceptanceCutLower = cutClower; fFMDCacceptanceCutUpper = cutCupper; }


    private:

        void                    PrintSetup();
        void                    CreateTHnCorrelations();

        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;

        Double_t                RangePhi(Double_t DPhi);
        Double_t                RangePhiFMD(Double_t DPhi);
        Double_t                GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius);
        void                    FillCorrelations(const Int_t spec);
        void                    FillCorrelationsMixed(const Int_t spec);
        Bool_t                  PrepareFMDTracks();

        Int_t                   IdentifyTrack(const AliAODTrack* track) const; // PID
        Bool_t                  HasTrackPIDTPC(const AliAODTrack* track) const; // is TPC PID OK for this track ?
        Bool_t                  HasTrackPIDTOF(const AliAODTrack* track) const; // is TOF PID OK for this track ?
        Bool_t                  AreEfficienciesLoaded();
        Double_t                GetEff(const Double_t dPt, const Int_t spec = 0, const Double_t dEta = 0.0);


        AliAnalysisTaskCorrForFlowFMD(const AliAnalysisTaskCorrForFlowFMD&); // not implemented
        AliAnalysisTaskCorrForFlowFMD& operator=(const AliAnalysisTaskCorrForFlowFMD&); // not implemented

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputListCharged;    //! output list
        TList*                  fInputListEfficiency;    //! input list
        TObjArray*              fTracksTrig[4]; //!
        TObjArray*              fTracksAss; //!
        AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
        AliPIDCombined*         fPIDCombined; //! AliPIDCombined container
        AliEventPoolManager*    fPoolMgr;  //!  event pool manager for Event Mixing
        //output histograms
        TH1D*                   fhEventCounter; //!
        TH1D*                   fhEventMultiplicity; //!
        TH3D*                   fhTrigTracks[4]; //!
        AliTHn*                 fhSE[4]; //!
        AliTHn*                 fhME[4]; //!
        TH1D*                   fhEfficiency[4]; //! not eta dependent
        TH1D*                   fhEfficiencyEta[4][4]; //! eta dependent (4 sectors)
        TH2D*                   fHistFMDeta; //! vs PVz

        //event and track selection
        AnaType                 fAnalType;
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fIsHMpp; // [kFALSE]
        Bool_t                  fDoPID; // [kFALSE]
        Bool_t                  fUseNch; // [kFALSE]
        Bool_t                  fUseEfficiency; // [kFALSE]
        Bool_t                  fEfficiencyEtaDependent; // [kFALSE]
        Bool_t                  fUseFMDcut; // [kTRUE]
        Bool_t                  fUseOppositeSidesOnly; // [kFALSE]
        UInt_t                  fFilterBit;
        Int_t                   fbSign;
        Int_t                   fNofTracks;
        Int_t                   fNchMin;
        Int_t                   fNchMax;
        Int_t                   fNOfSamples; //[1]
        Int_t                   fSampleIndex; //[0]
        Double_t                fPtMinTrig;
        Double_t                fPtMaxTrig;
        Double_t                fPtMinAss;
        Double_t                fPtMaxAss;
        Double_t                fFMDcutapar0; // [1.64755]
        Double_t                fFMDcutapar1; // [119.602]
        Double_t                fFMDcutcpar0; // [2.73426]
        Double_t                fFMDcutcpar1;  // [150.31]
        Double_t                fFMDAacceptanceCutLower; // [1.8]
        Double_t                fFMDAacceptanceCutUpper; // [4.8]
        Double_t                fFMDCacceptanceCutLower; // [1.8]
        Double_t                fFMDCacceptanceCutUpper; // [3.2]
        std::vector<Double_t>   fPtBinsTrigCharged;
        std::vector<Double_t>   fPtBinsAss;
        Double_t                fCentMin;
        Double_t                fCentMax;
        Double_t                fCentrality;
        Double_t                fAbsEtaMax;
        Double_t                fPVz;
        Double_t                fPVzCut; // [10.]
        Double_t                fPIDbayesPion; // [0.95]
        Double_t                fPIDbayesKaon; // [0.85]
        Double_t                fPIDbayesProton; // [0.85]
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

        ClassDef(AliAnalysisTaskCorrForFlowFMD, 4);
};

#endif
