/**************************************************************************
 *    Author:       Mikkel Petersen                                       *
 *    Framework for calculating di-hadron correlation                     *
 *    for extraction of v_n{2} coefficients of unidentified particles     *
 *    in TPC-FMD correlations while studying the eta decorrelation.       *
 *                                                                        *
 *    If used, modified, or distributed,                                  *
 *    please aknowledge the author of this code.                          *
 **************************************************************************/

#ifndef AliAnalysisTaskCorrForFlowEta_H
#define AliAnalysisTaskCorrForFlowEta_H

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


class AliAnalysisTaskCorrForFlowEta : public AliAnalysisTaskSE
{
    public:

                                AliAnalysisTaskCorrForFlowEta();
                                AliAnalysisTaskCorrForFlowEta(const char *name, Bool_t bUseEff, Bool_t bUseCalib);
        virtual                 ~AliAnalysisTaskCorrForFlowEta();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        enum                    AnaType { eTPCFMDA = 0, eTPCFMDC, eFMDAFMDC};
        enum                    ColSystem { sPP = 0, sPPb, sPbPb};

        // global flags
        void                    SetAnaType(AnaType type) { fAnalType = type; }
        void                    SetColSystem(ColSystem type) { fColSystem = type; }
        void                    SetIsMC(Bool_t mc = kTRUE, Bool_t tpc = kTRUE, Bool_t fmd = kTRUE) { fIsMC = mc; fIsTPCgen = tpc; fIsFMDgen = fmd; }
        void                    SetIsHMpp(Bool_t hm = kTRUE) { fIsHMpp = hm; fTrigger = AliVEvent::kHighMultV0; }
        void                    SetUseOppositeSidesOnly(Bool_t sides = kTRUE) { fUseOppositeSidesOnly = sides; }
        void                    SetSystematicsFlag(TString flag) { fSystematicsFlag = flag; }
        void                    SetSkipCorrelations(Bool_t flag = kTRUE) { fSkipCorr = flag; }
//        void                    SetIsAniparticleCheck(Bool_t flag = kTRUE, Bool_t antip = kTRUE) { fIsAntiparticleCheck = flag; fDoAntiparticleOnly = antip; }
        void                    SetRejectSecondariesFromMC(Bool_t flag = kTRUE) { fRejectSecondariesFromMC = flag; }
        void                    SetVetoJetEvents(Bool_t flag = kTRUE) { fVetoJetEvents = flag; }
        void                    SetJetEventsLowPtCut(Double_t cut) { fJetParticleLowPt = cut; }

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

        //track selection V0s
        void                    SetNSigmaTPC(Double_t cut) { fSigmaTPC = cut; }
        void                    SetnTPCcrossedRows(Int_t cut) { fnTPCcrossedRows = cut; }

        // correlation related
        void                    SetPtRangeTrig(Double_t min, Double_t max) {fPtMinTrig = min; fPtMaxTrig = max; }
        void                    SetCentBinsForMixing(Int_t nofBins, std::vector<Double_t> bins) { fNCentBins = nofBins; fCentBins = bins; }
        void                    SetNofSamples(Int_t n) { fNOfSamples = n; }

        // FMD related
        void                    SetUseFMDcut(Bool_t cut = kTRUE) { fUseFMDcut = cut; }
        void                    SetFMDcutParameters(Double_t par0a, Double_t par1a, Double_t par0c, Double_t par1c) { fFMDcutapar0 = par0a; fFMDcutapar1 = par1a; fFMDcutcpar0 = par0c; fFMDcutcpar1 = par1c; }
        void                    SetFMDacceptanceCuts(Double_t cutAlower, Double_t cutAupper, Double_t cutClower, Double_t cutCupper) { fFMDAacceptanceCutLower = cutAlower; fFMDAacceptanceCutUpper = cutAupper; fFMDCacceptanceCutLower = cutClower; fFMDCacceptanceCutUpper = cutCupper; }
        void                    SetBoostAMPT(Bool_t flag = kTRUE){ fBoostAMPT = flag; }

        // Bin related, TPC, FMD, dEta
        void                    SetBinsTPC(Int_t nBins, std::vector<Double_t> bins) { fNBinsTPCeta = nBins; fBinsTPCeta = bins; }
        void                    SetBinsFMD(Int_t nBins, std::vector<Double_t> bins) { fNBinsFMDeta = nBins; fBinsFMDeta = bins; }
        void                    SetNBinsdEta(Int_t NBins) { fNBinsdEta = NBins; }

    private:

        void                    PrintSetup();
        void                    DebugPrint();
        void                    CreateTHnCorrelations();

        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;

        Double_t                RangePhi(Double_t DPhi);
        Double_t                RangePhiFMD(Double_t DPhi);
        void                    FillCorrelations();
        void                    FillCorrelationsMixed();
        Bool_t                  PrepareTPCTracks();
        Bool_t                  PrepareFMDTracks();
        Bool_t                  PrepareMCTracks();

        Bool_t                  AreEfficienciesLoaded();
        Double_t                GetEff(const Double_t dPt, const Double_t dEta = 0.0);
        Int_t                   GetEtaRegion(const Double_t dEta);
        TString                 ReturnPPperiod(const Int_t runNumber) const;
        Double_t                TransverseBoost(const AliMCParticle *track);

        AliAnalysisTaskCorrForFlowEta(const AliAnalysisTaskCorrForFlowEta&); // not implemented
        AliAnalysisTaskCorrForFlowEta& operator=(const AliAnalysisTaskCorrForFlowEta&); // not implemented

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputListCharged;    //! output list
        TList*                  fInputListEfficiency;    //! input list
        TObjArray*              fTracksTrig; //!
        TObjArray*              fTracksAss; //!
        AliEventPoolManager*    fPoolMgr;  //!  event pool manager for Event Mixing
        //output histograms
        TH1D*                   fhEventCounter; //!
        TH1D*                   fhEventMultiplicity; //!
        AliTHn*                 fhTrigTracks; //!
        AliTHn*                 fhSE; //!
        AliTHn*                 fhME; //!
        TH2D*                   fhEfficiency; //! not eta dependent
        TH2D*                   fhEfficiencyEta[1][8]; //! eta dependent (8 sectors)
        TH2D*                   fHistFMDeta; //! vs PVz
        TH1D*                   fhCentCalib; //!
        TH1D*                   fhPT; //!
        TH2D*                   fh2FMDvsV0[4]; //!

        //event and track selection
        AnaType                 fAnalType;
        ColSystem               fColSystem;
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fIsMC; // [kFALSE]
        Bool_t                  fIsTPCgen; // [kFALSE]
        Bool_t                  fIsFMDgen; // [kFALSE]
        Bool_t                  fIsHMpp; // [kFALSE]
        Bool_t                  fUseNch; // [kFALSE]
        Bool_t                  fUseEfficiency; // [kFALSE]
        Bool_t                  fUseFMDcut; // [kTRUE]
        Bool_t                  fUseOppositeSidesOnly; // [kFALSE]
        Bool_t                  fUseCentralityCalibration; // [kFALSE]
        Bool_t                  fSkipCorr; // [kFALSE]
        Bool_t                  fVetoJetEvents; // [kFALSE]
        Bool_t                  fRejectSecondariesFromMC; // [kFALSE]
        Bool_t                  fBoostAMPT; // [kFALSE] = boost to CMS in pPb collisions for the gen level of AMPT
        UInt_t                  fFilterBit;
        Int_t                   fbSign;
        Int_t                   fRunNumber; // previous run
        Int_t                   fNofTracks;
        Int_t                   fNofMinHighPtTracksForRejection;
        Int_t                   fNchMin;
        Int_t                   fNchMax;
        Int_t                   fNBinsTPCeta;
        Int_t                   fNBinsFMDeta;
        Int_t                   fNBinsdEta; //[25]
        Int_t                   fnTPCcrossedRows; // [70]
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
        std::vector<Double_t>   fBinsTPCeta;
        std::vector<Double_t>   fBinsFMDeta;
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
        Double_t                fSigmaTPC; // [3.0]
        Double_t                fJetParticleLowPt; // [5.]
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

        ClassDef(AliAnalysisTaskCorrForFlowEta, 17);
};

#endif
