/**************************************************************************
 *    Author:       Zuzana Moravcova
 *    Modified by:  Debojit Sarkar                                 *
 *    Framework for calculating di-hadron correlation                     *
 *    for extraction of v_n{2} coefficients of identified particles       *
 *    including primary identified particles (pi, K, p)                   *
 *    and reconstructed "V0" particles (K0s, Lambda)                      *
 *    using TPC-TPC and TPC-FMD correlations.                             *
 *                                                                        *
 *    If used, modified, or distributed,                                  *
 *    please aknowledge the author of this code.                          *
 **************************************************************************/

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
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliPicoTrack.h"
#include "TLorentzVector.h"



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
        enum                    PartSpecies { eCharged = 0, ePion, eKaon, eProton, eK0s, eLambda, ePhi };
        enum                    ColSystem { sPP = 0, sPPb, sPbPb};

        // global flags
        void                    SetAnaType(AnaType type) { fAnalType = type; }
        void                    SetColSystem(ColSystem type) { fColSystem = type; }
        void                    SetDoPID(Bool_t pid = kTRUE) { fDoPID = pid; }
        void                    SetDoV0(Bool_t v0 = kTRUE) { fDoV0 = v0; }
	void                    SetDoPHI(Bool_t Phi = kTRUE) { fDoPHI = Phi; }
	void                    SetPHIkinematics(Bool_t phishift = kFALSE, Bool_t rapshift = kFALSE) { fshiftphi_PHI = phishift; fshiftrap_PHI = rapshift; }

        void                    SetIsMC(Bool_t mc = kTRUE, Bool_t tpc = kTRUE, Bool_t fmd = kTRUE) { fIsMC = mc; fIsTPCgen = tpc; fIsFMDgen = fmd; }
        void                    SetIsHMpp(Bool_t hm = kTRUE) { fIsHMpp = hm; }
        void                    SetUseOppositeSidesOnly(Bool_t sides = kTRUE) { fUseOppositeSidesOnly = sides; }
        void                    SetSystematicsFlag(TString flag) { fSystematicsFlag = flag; }
        void                    SetSkipCorrelations(Bool_t flag = kTRUE) { fSkipCorr = flag; }
        void                    SetIsAniparticleCheck(Bool_t flag = kTRUE, Bool_t antip = kTRUE) { fIsAntiparticleCheck = flag; fDoAntiparticleOnly = antip; }
        void                    SetRejectSecondariesFromMC(Bool_t flag = kTRUE) { fRejectSecondariesFromMC = flag; }
        void                    SetVetoJetEvents(Bool_t flag, Double_t cut, Double_t selectionval, Bool_t selectjetsinTPC) { fVetoJetEvents = flag; fJetParticleLowPt = cut; fJetvetoselectionval = selectionval; fselectjetsinTPC = selectjetsinTPC; }
	void                    SetParticlemassbias(Bool_t massbias = kFALSE, Bool_t massbias_Proton=kFALSE, Bool_t massbias_Lambda=kFALSE, Bool_t massbias_Phi=kFALSE) { fParticlemass_bias_corr = massbias; fcheckmassbias_Proton = massbias_Proton; fcheckmassbias_Lambda = massbias_Lambda; fcheckmassbias_Phi = massbias_Phi; }


        // event selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetPVZcut(Double_t cut) { fPVzCut = cut; }
        void                    SetUseNchRange(Bool_t range, Int_t min, Int_t max) { fUseNch = range; fNchMin = min; fNchMax = max; }
	void                    SetUseNchfor_eventmixing(Bool_t range) { fUseNchfor_eventmixing = range; }
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
        void                    SetV0dcaK0ToPVcut(Double_t cut) { fV0dcaK0ToPV = cut; }
        void                    SetV0dcaLambdaToPVcut(Double_t pos, Double_t neg) { fV0dcaPosLambdaToPV = pos; fV0dcaNegLambdaToPV = neg; }
        void                    SetV0dcaDaugters(Double_t cutK0, Double_t cutLam) { fV0dcaDaugtersK0 = cutK0; fV0dcaDaugtersLambda = cutLam; }
        void                    SetK0radius(Double_t min, Double_t max) { fK0radiusMin = min; fK0radiusMax = max; }
        void                    SetLambdaradius(Double_t min, Double_t max) { fLambdaradiusMin = min; fLambdaradiusMax = max; }
        void                    SetV0sCPAs(Double_t k0s, Double_t lambda) { fCutCPAK0s = k0s; fCutCPALambda = lambda; }
        void                    SetV0sTaus(Double_t k0s, Double_t lambda) { fCutTauK0s = k0s; fCutTauLambda = lambda; }
        void                    SetNSigmaTPC(Double_t cut) { fSigmaTPC = cut; }
	void                    SetNSigmaTPCTOF(Double_t cut) { fNSigmaTPCTOF = cut; }

        void                    SetnTPCcrossedRows(Int_t cut) { fnTPCcrossedRows = cut; }
        void                    SetMassRejWindowK0(Double_t cut) { fMassRejWindowK0 = cut; }
        void                    SetMassRejWindowLambda(Double_t cut) { fMassRejWindowLambda = cut; }
	void                    SetK0MassRange(Double_t min, Double_t max) { fMinK0Mass = min; fMaxK0Mass = max; }
	void                    SetLambdaMassRange(Double_t min, Double_t max) { fMinLambdaMass = min; fMaxLambdaMass = max; }
	void                    SetPhiMassRange(Double_t min, Double_t max) { fMinPhiMass = min; fMaxPhiMass = max; }


        // correlation related
        void                    SetPtRangeTrig(Double_t min, Double_t max) {fPtMinTrig = min; fPtMaxTrig = max; }
        void                    SetPtRangeAss(Double_t min, Double_t max) {fPtMinAss = min; fPtMaxAss = max; }
        void                    SetPtBins(std::vector<Double_t> bins) { fPtBinsTrigCharged = bins; }
        void                    SetPhiStarCur(Double_t phiStar) {fMergingCut = phiStar; }
        void                    SetCentBinsForMixing(Int_t nofBins, std::vector<Double_t> bins) { fNCentBins = nofBins; fCentBins = bins; }
        void                    SetNofSamples(Int_t n) { fNOfSamples = n; }
        void                    SetNofMbins(Int_t n) { fNbinsMinv = n; }
	void                    SetPvzBins(Int_t nofBins, std::vector<Double_t> bins)  {fNzVtxBins = nofBins; fzVtxBins = bins;}


        // FMD related
        void                    SetUseFMDcut(Bool_t cut = kTRUE) { fUseFMDcut = cut; }
        void                    SetFMDcutParameters(Double_t par0a, Double_t par1a, Double_t par0c, Double_t par1c) { fFMDcutapar0 = par0a; fFMDcutapar1 = par1a; fFMDcutcpar0 = par0c; fFMDcutcpar1 = par1c; }
        void                    SetFMDacceptanceCuts(Double_t cutAlower, Double_t cutAupper, Double_t cutClower, Double_t cutCupper) { fFMDAacceptanceCutLower = cutAlower; fFMDAacceptanceCutUpper = cutAupper; fFMDCacceptanceCutLower = cutClower; fFMDCacceptanceCutUpper = cutCupper; }
        void                    SetBoostAMPT(Bool_t flag = kTRUE){ fBoostAMPT = flag; }

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

        Int_t                   IdentifyTrack(const AliAODTrack* track); // PID
        void                    PrepareV0(); // V0
	void                    PreparePhi(); // Phi
        Bool_t                  IsV0(const AliAODv0* v0) const; // V0s selection
        Bool_t                  IsK0s(const AliAODv0* v0) const;
        Bool_t                  IsLambda(const AliAODv0* v0);
        Double_t                ProperLifetime(const AliAODv0* v0, const Double_t massPDG) const;
        Bool_t                  HasTrackPIDTPC(const AliAODTrack* track) const; // is TPC PID OK for this track ?
        Bool_t                  HasTrackPIDTOF(const AliAODTrack* track) const; // is TOF PID OK for this track ?
        Bool_t                  AreEfficienciesLoaded();
        Double_t                GetEff(const Double_t dPt, const Int_t spec = 0, const Double_t dEta = 0.0);
        Int_t                   GetEtaRegion(const Double_t dEta);
        TString                 ReturnPPperiod(const Int_t runNumber) const;
        Double_t                TransverseBoost(const AliMCParticle *track);

        AliAnalysisTaskCorrForFlowFMD(const AliAnalysisTaskCorrForFlowFMD&); // not implemented
        AliAnalysisTaskCorrForFlowFMD& operator=(const AliAnalysisTaskCorrForFlowFMD&); // not implemented

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputListCharged;    //! output list
        TList*                  fInputListEfficiency;    //! input list
        TObjArray*              fTracksTrig[7]; //!
	TObjArray*              fTracksTrig_Kaon_Phi; //!
        TObjArray*              fTracksAss; //!
        AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
        AliPIDCombined*         fPIDCombined; //! AliPIDCombined container
        AliEventPoolManager*    fPoolMgr;  //!  event pool manager for Event Mixing
        //output histograms
        TH1D*                   fhEventCounter; //!
        TH1D*                   fhEventMultiplicity; //!
	TH1D*                   fhEventMultiplicity_jetveto; //!
	TH1D*                   fhEventMultiplicity_massbias; //!
        AliTHn*                 fhTrigTracks[7]; //!
        AliTHn*                 fhSE[7]; //!
        AliTHn*                 fhME[7]; //!
        TH2D*                   fhEfficiency[6]; //! not eta dependent
        TH2D*                   fhEfficiencyEta[6][8]; //! eta dependent (8 sectors)
        TH2D*                   fHistFMDeta; //! vs PVz
        TH1D*                   fhV0Counter[3]; //!
        TH1D*                   fhK0sphi; //!
	TH1D*                   fhLambdaphi; //!
	TH1D*                   fhPhiphi; //!	   		
        TH1D*                   fhCentCalib; //!
        TH1D*                   fhPT[7]; //!
	TH1D*                   fhPT_trig[7]; //!
        TH2D*                   fhPTvsMinv[3]; //!
	TH2D*                   fhPTvsMinv_Phi_LS; //!
        TH2D*                   fh2FMDvsV0[4]; //!

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
	Bool_t                  fDoPHI; // [kFALSE]
	Bool_t                  fshiftphi_PHI; // [kFALSE]
        Bool_t                  fshiftrap_PHI; // [kFALSE]
        Bool_t                  fUseNch; // [kFALSE]
	Bool_t                  fUseNchfor_eventmixing; // [kFALSE]
        Bool_t                  fUseEfficiency; // [kFALSE]
        Bool_t                  fUseFMDcut; // [kTRUE]
        Bool_t                  fUseOppositeSidesOnly; // [kFALSE]
        Bool_t                  fUseCentralityCalibration; // [kFALSE]
        Bool_t                  fSkipCorr; // [kFALSE]
        Bool_t                  fIsAntiparticleCheck; // [kFALSE]
        Bool_t                  fDoAntiparticleOnly; // [kFALSE] == positive particles only and lambdas
        Bool_t                  fVetoJetEvents; // [kFALSE]
	Bool_t                  fselectjetsinTPC; // [kFALSE]
        Bool_t                  fRejectSecondariesFromMC; // [kFALSE]
        Bool_t                  fBoostAMPT; // [kFALSE] = boost to CMS in pPb collisions for the gen level of AMPT
        UInt_t                  fFilterBit;
        Int_t                   fbSign;
        Int_t                   fRunNumber; // previous run
        Double_t                fNofTracks;
        Int_t                   fNofMinHighPtTracksForRejection;
        Int_t                   fNchMin;
        Int_t                   fNchMax;
        Int_t                   fNbinsMinv; // [60]
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
        Double_t                fV0dcaK0ToPV; // [0.06]
        Double_t                fV0dcaNegLambdaToPV; // [0.25]
        Double_t                fV0dcaPosLambdaToPV; // [0.1]
        Double_t                fV0dcaDaugtersK0; // [1.]
        Double_t                fV0dcaDaugtersLambda; // [1.]
        Double_t                fK0radiusMin; // [0.5]
        Double_t                fK0radiusMax; // [200.]
        Double_t                fLambdaradiusMin; // [0.5]
        Double_t                fLambdaradiusMax; // [200.]
        Double_t                fCutCPAK0s; // [0.97]
        Double_t                fCutCPALambda; // [0.995]
        Double_t                fCutTauK0s; // [0.]
        Double_t                fCutTauLambda; // [0.]
        Double_t                fSigmaTPC; // [3.0]
	Double_t                fNSigmaTPCTOF; // [0.0]
	
        Double_t                fMassRejWindowK0; // [0.005]
        Double_t                fMassRejWindowLambda; // [0.01]
	Double_t                fJetvetoselectionval; //[0.5]
	
	Double_t                fMinK0Mass; // [0.44]
        Double_t                fMaxK0Mass; // [0.56]
        Double_t                fMinLambdaMass; // [1.08]
        Double_t                fMaxLambdaMass; // [1.15]
	Double_t                fMinPhiMass; // [0.99]
        Double_t                fMaxPhiMass; // [1.07]
	Bool_t                  fParticlemass_bias_corr;
	Bool_t                  fcheckmassbias_Proton;//(kFALSE),
        Bool_t                  fcheckmassbias_Lambda;//(kFALSE),
        Bool_t                  fcheckmassbias_Phi;//(kFALSE),
        Int_t                   fProtonSigcount;
        Int_t                   fLambdaSigcount;
        Int_t                   fPhiSigcount;


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

        ClassDef(AliAnalysisTaskCorrForFlowFMD, 18);
};

#endif
