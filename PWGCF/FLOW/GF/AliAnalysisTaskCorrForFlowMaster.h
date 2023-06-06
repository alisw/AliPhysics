#ifndef AliAnalysisTaskCorrForFlowMaster_H
#define AliAnalysisTaskCorrForFlowMaster_H

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
#include "AliStack.h"
#include "TVirtualMCStack.h"

class AliAnalysisTaskCorrForFlowMaster : public AliAnalysisTaskSE
{
    public:

                                AliAnalysisTaskCorrForFlowMaster();
                                AliAnalysisTaskCorrForFlowMaster(const char *name, Bool_t bUseEff, Bool_t bUseCalib);
        virtual                 ~AliAnalysisTaskCorrForFlowMaster();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        enum                    AnaType { eTPCFMDA = 0, eTPCFMDC, eFMDAFMDC, eTPCTPC };
        enum                    PartSpecies { eCharged = 0, ePion, eKaon, eProton, eK0s, eLambda };
        enum                    ColSystem { sPP = 0, sPPb, sPbPb};

        // global flags
        void                    SetAnaType(AnaType type) { fAnalType = type; }
        void                    SetColSystem(ColSystem type) { fColSystem = type; }
        void                    SetIsMC(Bool_t mc = kTRUE, Bool_t tpc = kTRUE) { fIsMC = mc; fIsTPCgen = tpc; }
        void                    SetUseEventBias(Bool_t eventbias = kTRUE) {fUseEventBias = eventbias;}
        void                    SetIsHMpp(Bool_t hm = kTRUE) { fIsHMpp = hm; }
        void                    SetUseOppositeSidesOnly(Bool_t sides = kTRUE) { fUseOppositeSidesOnly = sides; }
        void                    SetSystematicsFlag(TString flag) { fSystematicsFlag = flag; }
        void                    SetSkipCorrelations(Bool_t flag = kTRUE) { fSkipCorr = flag; }
        void                    SetRejectSecondariesFromMC(Bool_t flag = kTRUE) { fRejectSecondariesFromMC = flag; }
        void                    SetUsePhiStar(Bool_t flag) {fUsePhiStar = flag;}
        void                    SetUseEfficiency(Bool_t flag) {fUseEfficiency = flag;}
        void                    SetCreateQAPlots(Bool_t flag) {fCreateQAPlots = flag;}

        // event selection
        void                    SetNumEventBias(Int_t num) {fNumEventBias = num;}
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetPVZcut(Double_t cut) { fPVzCut = cut; }
        void                    SetUseNchRange(Bool_t range, Int_t min, Int_t max) { fUseNch = range; fNchMin = min; fNchMax = max; }
        void                    SetCentrality(TString cent, Double_t min = 0.0, Double_t max = 20.0) { fCentEstimator = cent; fCentMin = min; fCentMax = max; }
        void                    SetVetoJetEvents(Bool_t flag, Double_t cut, Double_t selectionval, Bool_t selectjetsinTPC) { fVetoJetEvents = flag; fJetParticleLowPt = cut; fJetvetoselectionval = selectionval; fselectjetsinTPC = selectjetsinTPC; }

        //track selection (charged + PID + som global )
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetAbsEta(Double_t etaAbs) { fAbsEtaMax = etaAbs; }
        void                    SetTPCclMincut(Double_t cut) { fTPCclMin = cut; }
        void                    SetutDCAz(Double_t cut) { fCutDCAz = cut; }
        void                    SetutDCAxy(Double_t cut) { fCutDCAxySigma = cut; }
        void                    SetCutTPCchi2pCl(Double_t cut) { fCutTPCchi2pCl = cut; }
        void                    SetUseLikeSign(Bool_t flag) {fUseLikeSign = flag;}
        void                    SetUseUnlikeSign(Bool_t flag) {fUseUnlikeSign = flag;}

        //track selection V0s
        void                    SetNSigmaTPC(Double_t cut) { fSigmaTPC = cut; }
        void                    SetnTPCcrossedRows(Int_t cut) { fnTPCcrossedRows = cut; }

        // correlation related
        void                    SetPtRangeTrig(Double_t min, Double_t max) {fPtMinTrig = min; fPtMaxTrig = max; }
        void                    SetPtRangeAss(Double_t min, Double_t max) {fPtMinAss = min; fPtMaxAss = max; }
        void                    SetPtBins(std::vector<Double_t> bins) { fPtBinsTrigCharged = bins; }
        void                    SetPtBinsAss(std::vector<Double_t> bins) { fPtBinsAssCharged = bins; }
        void                    SetPhiStarCur(Double_t phiStar) {fMergingCut = phiStar; }
        void                    SetCentBinsForMixing(Int_t nofBins, std::vector<Double_t> bins) { fNCentBins = nofBins; fCentBins = bins; }
        void                    SetNofSamples(Int_t n) { fNOfSamples = n; }
        void                    SetPtRefRange(Double_t min, Double_t max) {fPtRefMin=min; fPtRefMax=max; }

        //FMD
        void                    SetBoostAMPT(Bool_t flag = kTRUE){ fBoostAMPT = flag; }

    private:

        void                    PrintSetup();
        void                    CreateTHnCorrelations();

        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;

        Double_t                RangePhi(Double_t DPhi);
        Double_t                GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius);
        Bool_t                  CheckDPhiStar(Double_t dEta, Double_t trigPhi, Double_t trigPt, Double_t trigCharge, Double_t assPhi, Double_t assPt, Double_t assCharge);
        void                    FillCorrelations();
        void                    FillCorrelationsMixed();
        Bool_t                  PrepareTPCTracks();
        Bool_t                  PrepareMCTracks();

        Int_t                   IdentifyTrack(const AliAODTrack* track) const; // PID
        Bool_t                  AreEfficienciesLoaded();
        Double_t                GetEff(const Double_t dPt, const Int_t spec = 0, const Double_t dEta = 0.0);
        Int_t                   GetEtaRegion(const Double_t dEta);
        TString                 ReturnPPperiod(const Int_t runNumber) const;
        Double_t                TransverseBoost(const AliMCParticle *track);

        AliAnalysisTaskCorrForFlowMaster(const AliAnalysisTaskCorrForFlowMaster&); // not implemented
        AliAnalysisTaskCorrForFlowMaster& operator=(const AliAnalysisTaskCorrForFlowMaster&); // not implemented

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
        AliTHn*                 fhSEref; //!
        AliTHn*                 fhMEref; //!
        TH2D*                   fhEfficiency[1]; //! not eta dependent
        TH2D*                   fhEfficiencyEta[8]; //! eta dependent (8 sectors)
        TH1D*                   fhCentCalib; //!
        TH1D*                   fhPT; //!
        TH1D*                   fhPhi; //!
        TH1D*                   fhEta; //!
        TH1D*                   fhPVz; //!


        //event and track selection
        AnaType                 fAnalType;
        ColSystem               fColSystem;
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fIsMC; // [kFALSE]
        Bool_t                  fIsTPCgen; // [kFALSE]
        Bool_t                  fIsHMpp; // [kFALSE]
        Bool_t                  fUseEventBias; // [kFALSE]
        Bool_t                  fUseNch; // [kFALSE]
        Bool_t                  fUseEfficiency; // [kFALSE]
        Bool_t                  fUseOppositeSidesOnly; // [kFALSE]
        Bool_t                  fUseCentralityCalibration; // [kFALSE]
        Bool_t                  fSkipCorr; // [kFALSE]
        Bool_t                  fVetoJetEvents; // [kFALSE]
        Bool_t                  fRejectSecondariesFromMC; // [kFALSE]
        Bool_t                  fBoostAMPT; // [kFALSE] = boost to CMS in pPb collisions for the gen level of AMPT
        Bool_t                  fUsePhiStar; // [kFALSE]
        Bool_t                  fCreateQAPlots; //[kFALSE]
        Bool_t                  fUseLikeSign; // [kFALSE]
        Bool_t                  fUseUnlikeSign; // [kFALSE]
        Bool_t                  fselectjetsinTPC; // [kFALSE]
        UInt_t                  fFilterBit;
        Int_t                   fbSign;
        Int_t                   fRunNumber; // previous run
        Int_t                   fNofTracks;
        Int_t                   fNofTrackGlobal; // [0]
        Int_t                   fNofEventGlobal; // [0]
        Int_t                   fNofMinHighPtTracksForRejection;
        Int_t                   fNchMin;
        Int_t                   fNchMax;
        Int_t                   fNumEventBias; // [2]
        Int_t                   fnTPCcrossedRows; // [70]
        Double_t                fNOfSamples; //[1]
        Double_t                fSampleIndex; //[0]
        Double_t                fPtMinTrig;
        Double_t                fPtMaxTrig;
        Double_t                fPtMinAss;
        Double_t                fPtMaxAss;
        Double_t                fPtRefMin; // [0.2]
        Double_t                fPtRefMax; // [3.0]
        Double_t                fJetvetoselectionval; //[0.5]
        std::vector<Double_t>   fPtBinsTrigCharged;
        std::vector<Double_t>   fPtBinsAssCharged;
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

        ClassDef(AliAnalysisTaskCorrForFlowMaster, 17);
};

#endif