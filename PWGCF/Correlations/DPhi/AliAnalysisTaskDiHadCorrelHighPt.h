/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDiHadCorrelHighPt_H
#define AliAnalysisTaskDiHadCorrelHighPt_H

#include "AliAnalysisTaskSE.h"
#include "AliEventPoolManager.h"
#include "AliEventCuts.h"
#include "AliPIDCombined.h"

class AliPIDResponse;
class AliEventPoolManager;
class THnSparse;
class AliAODv0;
class AliAODVertex;
class AliESDEvent;
class AliESDtrack;
class AliESDcascade;
class AliMCParticle;
class THnSparse;

class AliAnalysisTaskDiHadCorrelHighPt : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskDiHadCorrelHighPt();
                                AliAnalysisTaskDiHadCorrelHighPt(const char *name, Bool_t analysisMC, Bool_t useeff);
        virtual                 ~AliAnalysisTaskDiHadCorrelHighPt();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        Bool_t                  IsMyGoodPrimaryTrack(const AliAODTrack* aodtrack);
        Bool_t                  IsMyGoodPrimaryTrackESD(const AliESDtrack* esdtrack);
        Bool_t                  IsMyGoodV0AngleLambda(const AliAODv0 *t, AliAODVertex *pv);
        Bool_t                  IsMyGoodV0AngleLambdaESD(const AliESDv0 *t);
        Bool_t                  IsMyGoodV0AngleK0(const AliAODv0 *t, AliAODVertex *pv);
        Bool_t                  IsMyGoodV0AngleK0ESD(const AliESDv0 *t);
        Bool_t                  IsMyGoodV0RapidityK0(const AliAODv0 *t);
        Bool_t                  IsMyGoodV0RapidityK0ESD(const AliESDv0 *t);
        Bool_t                  IsMyGoodV0RapidityLambda(const AliAODv0 *t);
        Bool_t                  IsMyGoodV0RapidityLambdaESD(const AliESDv0 *t);
        Bool_t                  IsMyGoodLifeTimeK0(const AliAODv0 *V0);
        Bool_t                  IsMyGoodLifeTimeLambda(const AliAODv0 *V0);
        Bool_t                  IsMyGoodLifeTimeESD(const AliESDv0 *V0, Double_t masshypotesis);
        Bool_t                  IsMyGoodLifeTimeAntiLambda(const AliAODv0 *V0);
        Bool_t                  IsMyGoodPID(const AliAODTrack *TrackPos, const AliAODTrack *TrackNeg);
        Bool_t                  IsMyGoodDaughterTrack(const AliAODTrack *t)  ;
        Bool_t                  IsMyGoodDaughterTrackESD(const AliESDtrack *t, const AliESDv0 *v0)  ;
        Bool_t                  IsMyGoodV0(const AliAODv0 *v0,const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg, Int_t oSta);
        Bool_t                  IsMyGoodV0ESD(const AliESDv0 *v0,const AliESDtrack* myTrackPos, const AliESDtrack* myTrackNeg, Int_t oSta);
        Bool_t                  IsMyGoodV0Topology(const AliAODv0 *v0, Bool_t K0s, Bool_t Lambda);
        Bool_t                  IsMyGoodV0TopologyESD(const AliESDv0 *v0, const AliESDtrack* myTrackPos, const AliESDtrack* myTrackNeg, Bool_t K0s, Bool_t Lambda);

        Int_t                   GetOStatus() { return fOStatus; }
        void                    SetOStatus(Int_t stat) {  fOStatus=stat; }
        void                    SetCutsCrosscheck(Bool_t stat) {  fCutsCrosscheck=stat; }

        void                    SetPtTrigMin(Double_t var) {fPtTrigMin=var;}
        void                    SetPtAsocMin(Double_t var) {fPtAsocMin=var;}
        void                    Corelations(TObjArray *triggers, TObjArray *associated,Bool_t hh,Bool_t V0h,Float_t perc,Bool_t hV0);
        void                    CorelationsMixing(TObjArray *triggers, TObjArray *bgTracks,Float_t perc);
        void                    TopologCuts(Double_t pttrig,Double_t mass,Double_t dcaNeg, Double_t dcaPos,Double_t dcaDau, Double_t V0rad, Double_t cosPA,Double_t lifetime,Double_t massSell,Double_t triggType, Double_t status);
        void                    FillMC(const AliVParticle *V0,Int_t pdgV0,Int_t pdgDau1, Int_t pdgDau2,Int_t triggerType, Double_t mass, Double_t posTrackProp[6],Double_t negTrackProp[6]);
        void                    CorelationsMixinghV0(TObjArray *bgTracks, TObjArray *assocArray, Float_t perc);

        void                    SetCosPAK0(Float_t cosPAK0) { fCosPointAngleK0 = cosPAK0; }
        void                    SetCosPALam(Float_t cosPALam) { fCosPointAngleLam = cosPALam; }
        void                    SetDCAV0Daughters(Float_t dcaK0, Float_t dcaLam) { fDCAV0DaughtersK0 = dcaK0; fDCAV0DaughtersLambda = dcaLam; }
        void                    SetDCADaughterK0(Float_t dcapos) { fDCADaughterK0 = dcapos; }
        void                    SetDCADaughterLambda(Float_t dca_pos, Float_t dca_neg) { fDCAposDaughterLamda = dca_pos; fDCAnegDaughterLamda = dca_neg; }
        void                    SetNTPCcrossedRows(Int_t Ncr) { fnumOfTPCcrossedRows = Ncr; }
        void                    SetEffAnalysis(Bool_t eff) { fEfficiency = eff; }
        void                    SetPurityCheckAnalysis(Bool_t pur) { fPurityCheck = pur; }
        void                    SetCorrelationsAnalysis(Bool_t corr) { fCorrelations = corr; }
        void                    SetMixedEvents(Int_t nmix) { fMixedEvents = nmix; }
        void                    SetNofPtBinsTrigger(Int_t nbins) { fNumberOfPtBinsTrigger = nbins; }
        void                    SetNofPtBinsAssoc(Int_t nbins) { fNumberOfPtBinsAssoc = nbins; }
        void                    SetV0RadiusCut(Double_t cutK0, Double_t cutLam) { fK0sRadius = cutK0; fLambdaRadius = cutLam; }
        void                    SetSigmaCut(Double_t sigma) { fSigmaCut = sigma; }
        void                    SetEtaCut(Double_t eta) { fEtaCut = eta; }
        void                    SetRapidityCut(Double_t rap) { fRapidityCut = rap; }
        void                    SetLifeTimeCutLambda(Double_t lifetime) { fLifeTimeLam = lifetime; }
        void                    SetLifeTimeCutK0(Double_t lifetime) { fLifeTimeK0 = lifetime; }
        void                    SetMassRejectCutK0(Double_t rej) { fMassRejectCutK0 = rej; }
        void                    SetLifeTimeCutLamda(Double_t rej) { fMassRejectCutLam = rej; }
        void                    SetEventPileUpCut(Bool_t cut) { fRejectEventPileUp = cut; }
        void                    SetTrackPileUpCut(Bool_t cut) { fRejectTrackPileUp = cut; }
        void                    SetTrackPileUpTOFCut(Bool_t cut) { fRejectTOF = cut; }
        void                    SetV0PileUpCut(Bool_t cut) { fRejectV0PileUp = cut; }
        void                    SetMultiplicityEstimator(TString est) { fMultEstimator = est; }
        void                    SetCorrelationsGen(Bool_t correl) { fCorrelationsGen = correl; }
        void                    SetCorrelationsV0h(Bool_t correl) { fV0hCorr = correl; }
        void                    SetCorrelationshh(Bool_t correl) { fhhCorr = correl; }
        void                    SetCorrelationshV0(Bool_t correl) { fhV0Corr = correl; }
        void                    SetCorrelationshPion(Bool_t correl) { fhPionCorr = correl; }
        void                    SetFilterBit (Int_t filter) { fFilterBit = filter; }
        void                    SetRemoveHadronsFromV0 (Bool_t rem) { fRemoveHadrFromV0 = rem; }
        void                    SetRemoveLamhFromCascade (Bool_t rem) { fRemoveLamhFromCascade = rem; }
        void                    SetAcceptLambdasFromCascades (Bool_t accept) { fAacceptLambdasFromCasscade = accept; }
        void                    SetAcceptPurePrimHadrons (Bool_t accept) { fPurePrimHadrons = accept; }
        void                    SetAcceptPureV0 (Bool_t accept) { fPureV0 = accept; }
        void                    SetNumberOfDeltaPhiBins(Int_t bins) { fNumberOfDeltaPhiBins = bins; }
        void                    SetNumberOfDeltaEtaBins(Int_t bins) { fNumberOfDeltaEtaBins = bins; }
        void                    SetNumberOfEtaBins(Int_t bins) { fNumberOfEtaBins = bins; }
        void                    SetMixing(Bool_t mix) { fMixing = mix; }
        void                    SetMixingGen (Bool_t mix) { fMixingGen = mix; }
        void                    SetNumberOfVzBins(Int_t nbins) { fNumOfVzBins = nbins; }
        AliAODTrack *           SetAliAODTrack(Double_t theta,Double_t phi, Double_t pt , Short_t charge);
        void                    SetPVCut(Int_t cut) { fPrimaryVertexCut = cut; }
        void                    SetESDtrackCuts(AliESDtrackCuts *const trcuts) { fESDTrackCuts = trcuts; }
        void                    SetAnalysisAOD(Bool_t aodAnalysis) { fAnalysisAOD = aodAnalysis; }
        void                    SetTestPions(Bool_t pions) { fTestPions = pions; }
        void                    SetAnalyseFeedDown(Bool_t fd) { fAnalyseFeedDown = fd; }
        void                    CorrelationsXi(TObjArray *triggers,TObjArray *assoc, Double_t perc,Bool_t Xih);
        void                    SetMaximalPtAssoc(Double_t max) { fPtAssocMax = max; }
        void                    SetMaximalPtTrigg(Double_t max) { fPtTrigMax = max; }
        Bool_t                  IsNotOOBPileUp(AliESDtrack* neg, AliESDtrack* pos, AliESDtrack* bach);
        Bool_t                  IsMyGoodXiCandidate(AliESDcascade * cas, Double_t * par);
        Bool_t                  IsMyGoodXiDecayRadius(AliESDcascade * cas);
        Bool_t                  IsMyGoodDaughterV0DecayRadius(AliESDcascade * cas);
        Bool_t                  IsMyGoodXiProperLifetime(AliESDcascade * cas);
        Bool_t                  IsMyGoodXiDaughterTracks(AliESDcascade * cas,Int_t casCharge, Int_t ids[3], Double_t massXi);
        Bool_t                  IsMyGoodXiDaughterTrackPID(AliESDtrack * tr, Int_t casCharge, Int_t trCharge);
        Bool_t                  IsMyGoodBachelorTrack(AliESDtrack * bach);
        Bool_t                  IsMyGoodXiDaughterV0Mass(AliESDcascade * cas, Int_t casChar);
        Int_t                   IsGoodMCV0(AliMCParticle *pos, AliMCParticle * neg);
        Int_t                   IsGoodMCCascade(AliMCParticle *v0, AliMCParticle * bach);
        void                    SetMaxDCAToVertexZ(Float_t max) { fMaxDCAToVertexZ = max; }
        void                    SetMaxDCAToVertexXY(Float_t max) { fMaxDCAToVertexXY = max; }
        void                    SetMinNCrossedRowsTPCprimtracks(Float_t min) { fMinNCrossedRowsTPCprimtracks = min; }
        void                    SetNPhiBins(Int_t nbins) { fNumberPhiBins = nbins; }
        void                    SetEffCorrMixing(Bool_t mix) { fMixCorrect = mix; }
        void                    SetMinBiasTrigger(Bool_t tr) { fminBias = tr; }
        void                    SetHighMultTriggerV0(Bool_t tr) { fhighMult = tr; }
        void                    SetHighMultTriggerSPD(Bool_t tr) { fhighMultSPD= tr; }
        void                    SetPercentileRange(Double_t min, Double_t max) { fPercentileMin = min; fPercetileMax = max; }
        void                    SetEventQAPlots(Bool_t pl) { fEventCutsQAPlots = pl; }
        void                    SetppAnalysis(Bool_t pp) { fpp = pp; }
        void                    SetNMultiplicityBins(Int_t nBins) { fNMultiplicityBins = nBins; }
        void                    SetCascadePointingAgle(Double_t pa) { fCascadePA = pa; }
        void                    SetCascadeTransverseDecayRadius(Double_t dr) { fCascadeTransverseDecayRadius = dr; }
        void                    SetDaughterV0DecayRadius(Double_t dr) { fDaughterV0DecayRadius = dr; }
        void                    SetDaughterV0cosPA(Double_t pa) { fDaughterV0cosPA = pa; }
        void                    SetDCADaughterV0PV(Double_t dca) { fDCADaughterV0PV = dca; }
        void                    SetDCADaughterV0tracks(Double_t dca) { fDCADaughterV0tracks = dca;}
        void                    SetDCACascadeDaughters(Double_t dca) { fDCACascadeDaughters = dca; }
        void                    SetDaugterV0InvariantMassRange(Double_t m) { fDaugterV0InvariantMassRange = m; }
        void                    SetCascadeProperLifetime(Double_t pl) { fCascadeProperLifetime = pl; }
        void                    SetNSigmasTPCCascadedaughterTracks(Double_t nsig) { fNSigmasTPCCascadedaughterTracks = nsig; }
        void                    SetDCAbachelorPV(Double_t dca) { fDCAbachelorPV = dca; }
        void                    SetDCAbaryonDaughterPV(Double_t dca) { fDCAbaryonDaughterPV = dca; }
        void                    SetDCAmesonDaughterPV(Double_t dca) { fDCAmesonDaughterPV = dca; }
        Int_t                   GetOriginalPartonPDG(const AliMCParticle* mcTrack);
        void                    SetOnTheFlyMCTrain(Bool_t tr) { fonTheFlyMC = tr; }
        void                    SetTPCrowsRindableRatio(Double_t ratio) { fTPCrowsRindableRatio = ratio; }
        void                    SetMinimalTrackLenght(Int_t cut) { fTrackLength = cut; }
        void                    SetCollisionSystem(TString system) { fSystem = system; }
        void                    SetBayesPIDPioncut(Double_t pion) { fPIDbayesPion = pion; }
        Bool_t                  IsPionTrack(const AliAODTrack* track);
        Bool_t                  HasTrackPIDTPC(const AliAODTrack* track);
        Bool_t                  HasTrackPIDTOF(const AliAODTrack* track);
        void                    LoadEfficiencies();
        void                    SetAxisTitles(THnSparse *hist);
        Double_t                GetEff(Double_t pt, Double_t eta,Int_t candidate);
        Int_t                   GetPvzBin();

        AliEventCuts *           fAliEventCuts; //!

    private:
        AliAODEvent*            fAOD;           		//! input event
        AliESDEvent*            fESD;                   //! input event
        AliMCEvent*             fmcEvent;               //! input MC event
        AliPIDResponse*         fPIDResponse;           //!
        AliPIDCombined*         fPIDCombined; //!
        TList*                  fOutputList;    		//! output list
        TH3F*					          fHistLambdaMassPtCut;	//!
        TH3F*					          fHistK0MassPtCut;		//!
        TH3F*                   fHistAntiLambdaMassPtCut; //!
        THnSparse*              fHistKorelacie;             //!
        THnSparse*              fHistdPhidEtaMix;           //!
        TH1D*                   fHistV0Multiplicity;         //!
        TH2D*                   fHistMultVtxz;              //!
        THnSparse*              fHistMCPtAs;            //!
        THnSparse*              fHistRCPtAs;            //!
        THnSparse*              fHistNumberOfTriggers;  //!
        THnSparse*              fHistMCKorelacie;       //!
        THnSparse*              fHistMCMixingRec;       //!
        THnSparse*              fHistMCMixingGen;       //!

        Bool_t          			  fFillMixed;  // enable event mixing (default: ON)
        Int_t           			  fMixingTracks;      // size of track buffer for event mixing
        AliEventPoolManager*          fPoolMgr;         //! event pool manager
        AliEventPool*                 fPool; //!
        AliEventPool*                 fPoolMCGen; //!
        Bool_t                        fAnalysisMC; // enable MC study
        Int_t                         fOStatus; //
        Double_t                      fPtTrigMin; //
        Double_t                      fPtAsocMin; //
        Bool_t                        fCutsCrosscheck; //
        Int_t                         fMixedEvents; // number of minimum mixed events
        Double_t                      fK0sRadius; // K0 radius
        Double_t                      fLambdaRadius; // Lambda radius
        Double_t                      fSigmaCut; // TPC PID
        Double_t                      fEtaCut; // pseudorapidity cut for primary hadrons
        Double_t                      fRapidityCut; // rapidity cut for V0
        Double_t                      fLifeTimeLam; // LifeTime cut for Lambdas and AnitLambdas
        Double_t                      fLifeTimeK0; // LifeTime cut for K0
        Double_t                      fMassRejectCutK0; // V0 competing rejection for K0
        Double_t                      fMassRejectCutLam;// V0 competing rejection for Lambdas and AntiLambdas
        Bool_t                        fRejectEventPileUp; // enable to use Pile-up cuts
        Bool_t                        fRejectTrackPileUp; // enable to use Bunch-Of Pile-up cuts for tracks

        THnSparse*              fHistKorelacieMCrec; //!
        THnSparse*              fHistNumberOfTriggersGen;  //!
        THnSparse*              fHistNumberOfTriggersRec;  //!

        THnSparse*              fHistRecV0; //!
        THnSparse*              fHistGenV0; //!

        TH3D*                   fHistRCPtTrigg;            //!
        TH1D*                   fHistSelection;            //!
        TH1F*                   fHistMultipPercentile;     //!
        THnSparse*              fHistTopolCut;             //!
        THnSparse*              fHistTopolCutMC;           //!
        THnSparse*              fHistPurityCheck;          //!

        Float_t                 fCosPointAngleK0; //
        Float_t                 fCosPointAngleLam; //
        Float_t                 fDCAV0DaughtersK0; //
        Float_t                 fDCAV0DaughtersLambda; //
        Float_t                 fDCADaughterK0; //
        Float_t                 fDCAposDaughterLamda; //
        Float_t                 fDCAnegDaughterLamda; //
        Int_t                   fnumOfTPCcrossedRows; //
        Bool_t                  fEfficiency; //
        Bool_t                  fPurityCheck;//
        Bool_t                  fCorrelations;//
        TH3F*                   fHistPtResolution;//!
        Int_t                   fNumberOfPtBinsTrigger; //
        Int_t                   fNumberOfPtBinsAssoc;//
        TH1F *                  fHistV0MultiplicityK0; //!
        TH1F *                  fHistV0Lam;//!
        TH1F *                  fHistMultiplicityALam;//!
        TH2F *                  fHitsNTracks;//!
        THnSparse *             fHistPhiEta;//!
        Bool_t                  fRejectTOF; // enable to use Bunch-Of Pile-up cuts for tracks using TOF
        Bool_t                  fRejectV0PileUp;// enable to use Bunch-Of Pile-up cuts for V0
        TString                 fMultEstimator; // enable to change the Multiplicity estimator
        Bool_t                  fCorrelationsGen; // enable to compute the correlation function from generated MC particles
        Bool_t                  fV0hCorr; // enable to run V0-h correlations separately
        Bool_t                  fhhCorr; // enable to run h-h correlations separately
        Bool_t                  fhV0Corr; // enable to run h-V0 correlations separately
        Bool_t                  fhPionCorr; //enable to run h-Pion correlations separately
        Int_t                   fFilterBit; // enable to vary filter bit for systematic studies
        Bool_t                  fRemoveLamhFromCascade; // enable to remove hh corelations, which are from the same V0
        Bool_t                  fRemoveHadrFromV0; // enable to remove Lamh corelations, which are from the same cascade
        Bool_t                  fAacceptLambdasFromCasscade; // accept Lambdas from cascades to efficiency calculation
        Bool_t                  fPurePrimHadrons; // anable to accept only pure primary hadrons for MC closure test
        Bool_t                  fPureV0; // anable to accept only pure good ID V0 for MC closure test
        Bool_t                  fAnalysisAOD;// anable to change between ESDs and AODs
        Int_t                   fNumberOfDeltaPhiBins; // Number of DeltaPhi Bins in correlation function and in mixing
        Int_t                   fNumberOfDeltaEtaBins; // Number of DeltaEta Bins in correlation functionand in mixing
        Int_t                   fNumberOfEtaBins; // Number of Eta bins for efficiency correction
        Int_t                   fNumberPhiBins; //Number of Phi bins for efficiency correction
        Bool_t                  fMixing; // enable mixing
        Bool_t                  fMixingGen; //
        Int_t                   fNumOfVzBins; // number of PV bins for mixing
        Int_t                   fPrimaryVertexCut; // PV position acceptance


        TH1D *                  fHistGenMultiplicity; //!
        TH2D *                  fHistTPCTracksVsClusters; //!
        TH2D *                  fHistVZeroPercentileTPCMult;//!

        AliESDtrackCuts *       fESDTrackCuts; //!
        Bool_t                  fTestPions; // for testing MC trains

        TH2D *                  fHistPosNegTracks; //!
        THnSparse *             fHistLambdaFeedDown; //!
        Bool_t                  fAnalyseFeedDown;
        Double_t                fPtAssocMax; // maximal pt of associated Particle
        Double_t                fPtTrigMax; // maximal pt of trigger Particle
        TH3F *                  fHistXiMinusMassPtCut; //!
        TH3F *                  fHistXiPlusMassPtCut; //!
        TH2D *                  fHistCasMC; //!

        Double_t                fMagneticField;
        Double_t                fPV[3];

        Float_t                 fMaxDCAToVertexZ; // DCA selection criterium for primary tracks in Z direction
        Float_t                 fMaxDCAToVertexXY; // DCA selection criterium for primary tracks in XY plane
        Float_t                 fMinNCrossedRowsTPCprimtracks; // TPC quality selection criterium for primary tracks

        TObjArray *             fmcTracksSel;  //!  generated associated particles
        TObjArray *             fmcGenTracksMixing; //! generated associated particles for Mixing
        TObjArray *             fmcTracksTrigSel; //! generated trigger charged hadrons
        TObjArray *             fmcTracksV0Sel; //! Generated V0 triggers
        TObjArray *             fmcV0AssocSel; //! Generated V0 assoc
        TObjArray *             fselectedMCassoc; //! all reconstructed associated particles, with reconstructed pt,phi,eta values - for raw correlation function
        TObjArray *             fselectedMCV0assoc; //! all reconstructed V0 as associated particles, with reconstructed pt,phi,eta values - for raw correlation function
        TObjArray *             fselectedMCtrig; //! all reconstructed trigger particles, with reconstructed pt,phi,eta values - for raw correlation function
        TObjArray *             fselectedMCV0Triggersrec; //! All reconstructed V0 candidates for triggers with reconstructed pt,phi,eta values - for raw correlation function
        TObjArray *             fselectedTracks; //! tracks for mixing
        TObjArray *             fselectedAssociatedTracks; //!
        TObjArray *             fselectedTriggerTracks; //!
        TObjArray *             fselectedV0Triggers; //!
        TObjArray *             fselectedV0Assoc; //!
        TH2D *                  fHistEffCorrectionHadron[9]; //!
        TH2D *                  fHistEffCorrectionK0[9]; //!
        TH2D *                  fHistEffCorrectionLam[9]; //!
        TH2D *                  fHistEffCorrectionAntiLam[9]; //!
        TH2D *                  fHistEffCorrectionNegXi[9]; //!
        TH2D *                  fHistEffCorrectionPosXi[9]; //!
        TH2D *                  fHistEffCorrectionNegPion[9]; //!
        TH2D *                  fHistEffCorrectionPosPion[9]; //!
        TH1F *                  fHistSecondaryCont; //!
        TH1F *                  fHistSecondaryContPosPion; //!
        TH1F *                  fHistSecondaryContNegPion; //!
        TList *                 fEffList; //!
        Bool_t                  fUseEff;
        Bool_t                  fMixCorrect; // enable efficiency correction for mixing
        Bool_t                  fminBias;// enable minimum bias trigger , deffault
        Bool_t                  fhighMult;// enable V0 high multiplicity trigger
        Bool_t                  fhighMultSPD;// enable SPD high multiplicity trigger
        TH3D *                  fHistMultVZEROTracklets; //!
        Double_t                fPercentileMin; // minimum muliplicity percentile to be analysed
        Double_t                fPercetileMax; // maximum muliplicity percentile to be analysed
        Bool_t                  fEventCutsQAPlots; // enable to save event QA plots
        Bool_t                  fpp; // flag pp / PbPb
        Int_t                   fNMultiplicityBins; // number of multiplicity bins
        Double_t                fCascadePA; //casscade cos(PA) cut
        Double_t                fCascadeTransverseDecayRadius; // cascade transverse radius cut
        Double_t                fDaughterV0DecayRadius; // cascade daughter V0 transverse decay radius cut
        Double_t                fDaughterV0cosPA; // cascade daughter V0 cos(PA) cut
        Double_t                fDCADaughterV0PV; // cascade DCA daughter V0 to PV cut
        Double_t                fDCADaughterV0tracks; // cascade DCA of daughter V0 tracks cut
        Double_t                fDCACascadeDaughters; // DCA of Cascade Daughters cut
        Double_t                fDaugterV0InvariantMassRange; // cascade daugter V0 invariant mass range cut
        Double_t                fCascadeProperLifetime; // cascade proper lifetime cut
        Double_t                fNSigmasTPCCascadedaughterTracks; // N Sigmas TPC of Cascade daughter Tracks cut
        Double_t                fDCAbachelorPV; // DCA of cacscade bachelor track to PV cut
        Double_t                fDCAbaryonDaughterPV; // DCA of cascade baryon daughter track to PV cut
        Double_t                fDCAmesonDaughterPV; // DCA of cascade meson daughter track to PV cut
        TH1D *                  fHistV0AmplitudeVsPVposition; //!
        Bool_t                  fonTheFlyMC; // on the fly MC train usage
        Double_t                fPercentile;// th V0M multiplicity percentile
        Double_t                fTPCrowsRindableRatio; // ratio of crossed rows in TPC and Finadable clusters
        Int_t                   fTrackLength;//V0 daughter track lenght
        TString                 fSystem; //
        Double_t                fPIDbayesPion; // [0.95]

        AliAnalysisTaskDiHadCorrelHighPt(const AliAnalysisTaskDiHadCorrelHighPt&); // not implemented
        AliAnalysisTaskDiHadCorrelHighPt& operator=(const AliAnalysisTaskDiHadCorrelHighPt&); // not implemented

        ClassDef(AliAnalysisTaskDiHadCorrelHighPt, 44);
};

class AliV0ChParticle : public AliVParticle
{
    public:
      AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate,Double_t mass)
        : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(-1),fIDh(-1),fIDpos(-1),fIDneg(-1),fIDbach(0),fMass(mass),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(-30)
      {
      }
      AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label)
        : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label),fIDh(-1),fIDpos(-1),fIDneg(-1),fIDbach(0),fMass(0),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(-30)
      {
      }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t iDh,Short_t charge)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh),fIDpos(-1),fIDneg(-1),fIDbach(0),fMass(0),fCharge(charge),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(-30)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t iDh,Int_t partontype)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh),fIDpos(-1),fIDneg(-1),fIDbach(0),fMass(0),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(partontype)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Short_t charge, Double_t pz, Double_t energ)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1),fIDpos(-1),fIDneg(-1),fIDbach(0),fMass(0),fCharge(charge), fPz(pz), fEnergie(energ),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(-30)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t iDh,Short_t charge, Double_t pz, Double_t energ)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh),fIDpos(-1),fIDneg(-1),fIDbach(0),fMass(0),fCharge(charge), fPz(pz), fEnergie(energ),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(-30)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t idpos, Int_t idneg,Double_t mass,Double_t phi1, Double_t pt1,Double_t eta1, Double_t char1, Double_t phi2,Double_t pt2,Double_t eta2, Double_t char2)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1),fIDpos(idpos), fIDneg(idneg), fIDbach(0),fMass(mass),fCharge(0),fPz(0),fEnergie(0),fPhi1(phi1),fPt1(pt1),fEta1(eta1),fChar1(char1),fPhi2(phi2),fPt2(pt2),fEta2(eta2),fChar2(char2),fOriginalParton(-30)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t idpos, Int_t idneg,Double_t mass)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1),fIDpos(idpos), fIDneg(idneg), fIDbach(0),fMass(mass),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(-30)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t idpos, Int_t idneg,Double_t mass, Int_t partontype)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1),fIDpos(idpos), fIDneg(idneg), fIDbach(0),fMass(mass),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(partontype)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt,Short_t candidate,Double_t mass,Int_t idpos, Int_t idneg,Int_t idbach)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(0), fIDh(-1),fIDpos(idpos), fIDneg(idneg), fIDbach(idbach),fMass(mass),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0),fOriginalParton(-30)
    {

    }
     virtual ~AliV0ChParticle() {}

      // kinematics
      virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Pz() const { return fPz; }
      virtual Double_t Pt() const { return fpT; }
      virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
      virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

      virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
      virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

      virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Phi()        const { return fPhi; }
      virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }


      virtual Double_t E()          const { return fEnergie; }
      virtual Double_t M()          const { return fMass; }

      virtual Double_t Eta()        const { return fEta; }
      virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }

      virtual Short_t Charge()      const { return fCharge; }
      virtual Int_t   GetLabel()    const { return fLabel; }
      // PID
      virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
      virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

      virtual Short_t WhichCandidate()      const { return fCandidate; }
    Int_t   GetIDCh()             const { return fIDh; }
    Int_t   GetIDPos()            const { return fIDpos; }
    Int_t   GetIDNeg()            const { return fIDneg; }
    Int_t   GetIDBach()            const { return fIDbach; }
    Short_t GetCharge1()          const { return fChar1; }
    Short_t GetCharge2()          const { return fChar2; }
    Double_t GetPhi1()            const { return fPhi1; }
    Double_t GetPhi2()            const { return fPhi2; }
    Double_t GetPt1()             const { return fPt1; }
    Double_t GetPt2()             const { return fPt2; }
    Double_t GetEta1()            const { return fEta1; }
    Double_t GetEta2()            const { return fEta2; }
    Int_t GetOrigalPartonType()   const { return fOriginalParton; }

      private:
      Double_t fEta;      // eta
      Double_t fPhi;      // phi
      Double_t fpT;       // pT
      Short_t fCandidate;   // V0 candidate: 1 - K0, 2 - Lambda, 3 - Antilambda, 4 - ChTrack
      Int_t fLabel;   // Label of MC particles
      Int_t fIDh;   // Label
      Int_t fIDpos;   // Label of possitive charged daughter
      Int_t fIDneg;   // Label of negative charged daughter
      Int_t fIDbach;   // Label of bachelor track
      Double_t fMass; // mass
      Short_t fCharge; // charge of the track
      Double_t fPz; //pZ
      Double_t fEnergie; //E
      Double_t fPhi1; // phi of first daughter
      Double_t fPt1; //pt of first daughter
      Double_t fEta1; // eta of first daughter
      Double_t fChar1; // charge of first daughter
      Double_t fPhi2; // phi of second daughter
      Double_t fPt2; //pt of second daughter
      Double_t fEta2; // eta of second daughter
      Double_t fChar2; // charge of second daughter
      Int_t fOriginalParton; //0- quark; 1-gluon, for MC truth information

      ClassDef( AliV0ChParticle, 8) // class required for correlatios calculation and event mixing
 };

#endif
