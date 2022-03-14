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
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliEventPoolManager.h"
#include "AliVEvent.h"
#include "AliTHn.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODv0.h"
#include "AliAODForwardMult.h"
#include "AliAODVZERO.h"
#include "AliPartSimpleForCorr.h"

class AliAnalysisTaskCorrForFlowFMD : public AliAnalysisTaskSE
{
    public:

                                AliAnalysisTaskCorrForFlowFMD();
                                AliAnalysisTaskCorrForFlowFMD(const char *name, Bool_t bUseEff, Bool_t bUseCalib);
        virtual                 ~AliAnalysisTaskCorrForFlowFMD();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        enum                    AnaType { eTPCFMDA = 0, eTPCFMDC, eFMDAFMDC, eTPCTPC };
        enum                    PartSpecies { eCharged = 0, ePion, eKaon, eProton, eK0s, eLambda };
        enum                    ColSystem { sPP = 0, sPPb};

        // global flags
        void                    SetAnaType(AnaType type) { fAnalType = type; }
        void                    SetColSystem(ColSystem type) { fColSystem = type; }
        void                    SetDoPID(Bool_t pid = kTRUE) { fDoPID = pid; }
        void                    SetDoV0(Bool_t v0 = kTRUE) { fDoV0 = v0; }
        void                    SetIsMC(Bool_t mc = kTRUE, Bool_t tpc = kTRUE, Bool_t fmd = kTRUE) { fIsMC = mc; fIsTPCgen = tpc; fIsFMDgen = fmd; }
        void                    SetIsHMpp(Bool_t hm = kTRUE) { fIsHMpp = hm; }
        void                    SetUseOppositeSidesOnly(Bool_t sides = kTRUE) { fUseOppositeSidesOnly = sides; }
        void                    SetSystematicsFlag(TString flag) { fSystematicsFlag = flag; }

        // event selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetPVZcut(Double_t cut) { fPVzCut = cut; }
        void                    SetUseNchRange(Bool_t range, Int_t min, Int_t max) { fUseNch = range; fNchMin = min; fNchMax = max; }
        void                    SetCentrality(TString cent, Double_t min = 0.0, Double_t max = 20.0) { fCentEstimator = cent; fCentMin = min; fCentMax = max; }

        //track selection (charged + PID + som global )
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetAbsEta(Double_t etaAbs) { fAbsEtaMax = etaAbs; }
        void                    SetTPCclMincut(Double_t cut) { fTPCclMin = cut; }
        void                    SetutDCAz(Double_t cut) { fCutDCAz = cut; }
        void                    SetutDCAxy(Double_t cut) { fCutDCAxySigma = cut; }
        void                    SetCutTPCchi2pCl(Double_t cut) { fCutTPCchi2pCl = cut; }
        void                    SetBayesPIDcut(Double_t pion, Double_t kaon, Double_t proton){ fPIDbayesPion = pion; fPIDbayesKaon = kaon; fPIDbayesProton = proton; }

        //track selection V0s
        void                    SetTPCclV0Ratio(Double_t cut) { fV0ratioClusters = cut; }
        void                    SetV0dcaToPVcut(Double_t cut) { fV0dcaToPV = cut; }
        void                    SetV0dcaDaugters(Double_t cut) { fV0dcaDaugters = cut; }
        void                    SetV0radius(Double_t min, Double_t max) { fV0radiusMin = min; fV0radiusMax = max; }
        void                    SetV0sCPAs(Double_t k0s, Double_t lambda) { fCutCPAK0s = k0s; fCutCPALambda = lambda; }
        void                    SetV0sTaus(Double_t k0s, Double_t lambda) { fCutTauK0s = k0s; fCutTauLambda = lambda; }

        // correlation related
        void                    SetPtRangeTrig(Double_t min, Double_t max) {fPtMinTrig = min; fPtMaxTrig = max; }
        void                    SetPtRangeAss(Double_t min, Double_t max) {fPtMinAss = min; fPtMaxAss = max; }
        void                    SetPtBins(std::vector<Double_t> bins) { fPtBinsTrigCharged = bins; }
        void                    SetPhiStarCur(Double_t phiStar) {fMergingCut = phiStar; }
        void                    SetCentBinsForMixing(Int_t nofBins, std::vector<Double_t> bins) { fNCentBins = nofBins; fCentBins = bins; }
        void                    SetNofSamples(Int_t n) { fNOfSamples = n; }
        void                    SetNofMbins(Int_t n) { fNbinsMinv = n; }

        // FMD related
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
        Bool_t                  PrepareTPCTracks();
        Bool_t                  PrepareFMDTracks();
        Bool_t                  PrepareMCTracks();

        Int_t                   IdentifyTrack(const AliAODTrack* track) const; // PID
        void                    PrepareV0(); // V0
        Bool_t                  IsV0(const AliAODv0* v0) const; // V0s selection
        Bool_t                  IsK0s(const AliAODv0* v0) const;
        Bool_t                  IsLambda(const AliAODv0* v0) const;
        Double_t                ProperLifetime(const AliAODv0* v0, const Double_t massPDG) const;
        Bool_t                  HasTrackPIDTPC(const AliAODTrack* track) const; // is TPC PID OK for this track ?
        Bool_t                  HasTrackPIDTOF(const AliAODTrack* track) const; // is TOF PID OK for this track ?
        Bool_t                  AreEfficienciesLoaded();
        Double_t                GetEff(const Double_t dPt, const Int_t spec = 0, const Double_t dEta = 0.0);
        Int_t                   GetEtaRegion(const Double_t dEta);
        TString                 ReturnPPperiod(const Int_t runNumber) const;


        AliAnalysisTaskCorrForFlowFMD(const AliAnalysisTaskCorrForFlowFMD&); // not implemented
        AliAnalysisTaskCorrForFlowFMD& operator=(const AliAnalysisTaskCorrForFlowFMD&); // not implemented

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputListCharged;    //! output list
        TList*                  fInputListEfficiency;    //! input list
        TObjArray*              fTracksTrig[6]; //!
        TObjArray*              fTracksAss; //!
        AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
        AliPIDCombined*         fPIDCombined; //! AliPIDCombined container
        AliEventPoolManager*    fPoolMgr;  //!  event pool manager for Event Mixing
        //output histograms
        TH1D*                   fhEventCounter; //!
        TH1D*                   fhEventMultiplicity; //!
        AliTHn*                 fhTrigTracks[6]; //!
        AliTHn*                 fhSE[6]; //!
        AliTHn*                 fhME[6]; //!
        TH2D*                   fhEfficiency[4]; //! not eta dependent
        TH2D*                   fhEfficiencyEta[4][8]; //! eta dependent (8 sectors)
        TH2D*                   fHistFMDeta; //! vs PVz
        TH1D*                   fhV0Counter[2]; //!
        TH1D*                   fhCentCalib; //!

        //event and track selection
        AnaType                 fAnalType;
        ColSystem               fColSystem;
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fIsMC; // [kFALSE]
        Bool_t                  fIsTPCgen; // [kFALSE]
        Bool_t                  fIsFMDgen; // [kFALSE]
        Bool_t                  fIsHMpp; // [kFALSE]
        Bool_t                  fDoPID; // [kFALSE]
        Bool_t                  fDoV0; // [kFALSE]
        Bool_t                  fUseNch; // [kFALSE]
        Bool_t                  fUseEfficiency; // [kFALSE]
        Bool_t                  fUseFMDcut; // [kTRUE]
        Bool_t                  fUseOppositeSidesOnly; // [kFALSE]
        Bool_t                  fUseCentralityCalibration; // [kFALSE]
        UInt_t                  fFilterBit;
        Int_t                   fbSign;
        Int_t                   fRunNumber; // previous run
        Int_t                   fNofTracks;
        Int_t                   fNchMin;
        Int_t                   fNchMax;
        Int_t                   fNbinsMinv; // [60]
        Double_t                fNOfSamples; //[1]
        Double_t                fSampleIndex; //[0]
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
        Double_t                fCentMin;
        Double_t                fCentMax;
        Double_t                fCentrality;
        Double_t                fAbsEtaMax;
        Double_t                fPVz;
        Double_t                fPVzCut; // [10.]
        Double_t                fTPCclMin; // [70.]
        Double_t                fCutDCAz; // [0.]
        Double_t                fCutDCAxySigma; // [0.]
        Double_t                fCutTPCchi2pCl; // [0.]
        Double_t                fPIDbayesPion; // [0.95]
        Double_t                fPIDbayesKaon; // [0.85]
        Double_t                fPIDbayesProton; // [0.85]
        Double_t                fV0ratioClusters; // [0.8]
        Double_t                fV0dcaToPV; // [0.06]
        Double_t                fV0dcaDaugters; // [1.]
        Double_t                fV0radiusMin; // [0.5]
        Double_t                fV0radiusMax; // [200.]
        Double_t                fCutCPAK0s; // [0.97]
        Double_t                fCutCPALambda; // [0.995]
        Double_t                fCutTauK0s; // [0.]
        Double_t                fCutTauLambda; // [0.]
        TString                 fCentEstimator; //"V0M"
        TString                 fSystematicsFlag; // ""
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

        ClassDef(AliAnalysisTaskCorrForFlowFMD, 9);
};

#endif
