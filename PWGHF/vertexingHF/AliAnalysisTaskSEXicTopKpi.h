#ifndef ALIANALYSISTASKSEXICTOPKPI_H
#define ALIANALYSISTASKSEXICTOPKPI_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

///*************************************************************************
/// \class Class AliAnalysisTaskSEXicTopKpi
/// \brief AliAnalysisTaskSE for Xic->pKpi candidates invariant mass histogram
/// and comparison to MC truth (kinematics stored in the AOD) and cut variables
/// distributions
/// \author Authors: A.Rossi, andrea.rossi@cern.ch
/// \author and M.Faggin, mattia.faggin@cern.ch
///*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <THnSparse.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsXictopKpi.h"
#include "AliNormalizationCounter.h"
#include "AliAODMCParticle.h"
#include "AliVertexerTracks.h"
#include "AliAODMCHeader.h"

class AliAODEvent;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliPIDResponse;
class AliAODVertex;
class AliAODRecoDecayHF3Prong;
class AliAnalysisVertexingHF;

class AliAnalysisTaskSEXicTopKpi : public AliAnalysisTaskSE
{
 public:
  enum ECandStatus {kGenLimAcc=1,kGenAccMother,kGenAcc,kReco=6,kRecoCuts,kRecoPID,kRecoLc=13,kRecoLcCuts,kRecoLcPID};

  AliAnalysisTaskSEXicTopKpi();
  AliAnalysisTaskSEXicTopKpi(const char *name,AliRDHFCutsD0toKpi* cuts);
  virtual ~AliAnalysisTaskSEXicTopKpi();


  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

/*   void CreateMCAcceptanceHistos(); */
/*   Bool_t CheckAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau); */
/*   void FillMCAcceptanceHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader); */

/*   void NormIPvar(AliAODEvent *aod, AliAODRecoDecayHF2Prong *part,TClonesArray *arrMC); */
/*   void SetArray(Int_t type=AliAnalysisTaskSEXicTopKpi::kD0){fArray=type;} */
/*   enum{kD0,kLS}; */
  
  
  AliESDtrack* SelectTrack(AliAODTrack *aodtr, Int_t &isSelProton,Int_t &isSelKaon, Int_t &isSelPion,Int_t &isSelSoftPion,AliESDtrackCuts *cutsProton, AliESDtrackCuts *cutsKaon, AliESDtrackCuts *cutsPion,AliESDtrackCuts *cutsSoftPion);

  void IsSelectedPID(AliAODTrack *track,Int_t &iSelPion,Int_t &iSelKaon,Int_t &iSelProton,const Int_t iSelPionCuts=1,const Int_t iSelKaonCuts=1,const Int_t iSelProtonCuts=1,Bool_t fillHistos=kFALSE);

  void SetReadMC(Bool_t readMC=kFALSE){fReadMC=readMC;}
  void SetAnalysisType(Int_t antype){fAnalysisType=antype;};
  void SetAODMismatchProtection(Int_t opt=0) {fAODProtection=opt;} 
  //void SetLcCuts(AliRDHFCutsLctopKpi *cuts){fCutsLc=cuts;}
  void SetXicCuts(AliRDHFCutsXictopKpi *cuts){fCutsXic=cuts;}
  Int_t CheckXicpKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab)const;
  void SetRecalcOwnPrimVtx(Bool_t recVtx){fRecalPrimVtx=recVtx;}
  Bool_t GetIsRecalcOwnPrimVtx(){return fRecalPrimVtx;}
  void SetSystem(Int_t sys){fSys=sys;}
  Int_t GetSystem(){return fSys;}
  void FillDist12and23(AliAODRecoDecayHF3Prong *pr,Double_t magfield);
  void SetUseLcTrackFilteringCut(Bool_t useLcTrackFilteringCut){fSetTrackCutLcFilteringPP=useLcTrackFilteringCut;}
  Int_t FlagCandidateWithVariousCuts(AliAODRecoDecayHF3Prong *pr,AliAODEvent *aod,Int_t massHypo);
  void SetMaxPtSPDkFirst(Bool_t applykfirst,Double_t minpt){
    fApplykFirst=applykfirst;
    fMaxPtTrackkFirst=minpt;
  }
  void SetFillTree(Int_t filltree){fFillTree=filltree;}
  void FillTree(AliAODRecoDecayHF3Prong *cand,Int_t massHypothesis,Float_t *varPointer,Int_t flagMC,AliAODEvent *aod,AliAODMCParticle* p,TClonesArray* array_MC, AliAODMCHeader *mcHeader);
  void SetMaxChi2Cut(Double_t maxchi2){fMaxVtxChi2Cut=maxchi2;}
  Double_t GetMaxChi2Cut(){return fMaxVtxChi2Cut;}
  Double_t CosThetaStar(Double_t mumVector[3],Double_t daughtVector[3],Double_t massMum,Double_t massDaught);
  void PrintCandidateVariables(AliAODRecoDecayHF3Prong *d,AliAODEvent *aod);  
  // downsampling for fTreeVar filling
  void SetDownsampling(Float_t pT_thr, Float_t down_lowpT, Float_t down_highpT){
    fpT_down = pT_thr;
    fLowpT_down = down_lowpT;
    fHighpT_down = down_highpT;
  }
  // pT limits for TTree filling
   void SetpTlimsTTreeFilling(Float_t min, Float_t max){
     fminpT_treeFill = min;
     fmaxpT_treeFill = max;
   }
  // require the calculation of dist12 and dist23
  void SetCalculate_dist12_dist23(Bool_t flag){ fCompute_dist12_dist23 = flag; }
  Short_t SetMapCutsResponse(Int_t massHypo_filtering, Int_t response_onlyCuts, Int_t response_onlyPID);

  // exporation of PID cuts with standard strategy
  void SetExplorePIDstd(Bool_t flag){ fExplore_PIDstdCuts=flag; }
  // dirty solution: flag to reduce the axes in the reco sparses ---> make the merging easier (mfaggin)
  void SetOnlyBayesPIDbin_recoSparse(Bool_t flag) {fOnlyBayesPIDbin=flag;}
  /// include the PID selection with Bayes approach only for proton
  void SetExplPID_BayesOnlyProt(Bool_t flag, Bool_t rejectStdPIDcases=kFALSE){
    if(!fExplore_PIDstdCuts)  SetExplorePIDstd(flag);
    fExplPID_BayesOnlyProt = flag;
    fNoStdPIDcases = rejectStdPIDcases;
  }

  void SetLcMassWindowForSigmaC(Double_t massrange){fLcMassWindowForSigmaC=massrange;}
  void SetSigmaCDeltaMassWindow(Double_t maxDeltaM){fSigmaCDeltaMassWindow=maxDeltaM;}
  void SetOnTheFlyLcCandidatesForSigmaC(Bool_t onthefly){fSigmaCfromLcOnTheFly=onthefly;}
  void SetFillOnlyTrackSparse(Bool_t fillonlysparse){fCheckOnlyTrackEfficiency=fillonlysparse;}
  void SetIsCdeuteronAnalysis(Bool_t iscd){fIsCdeuteronAnalysis=iscd;}
  void SetIsXicUpgradeAnalysis(Bool_t flag){fIsXicUpgradeAnalysis=flag;}
  void SetKeepOnlySigXicUpgradeAnalysis(Bool_t flag){fIsKeepOnlySigXicUpgradeAnalysis=flag;}
  void SetKeepOnlyBkgXicUpgradeAnalysis(Bool_t flag){fIsKeepOnlyBkgXicUpgradeAnalysis=flag;}
  void SetRejFactorBkgUpgrade(Double_t value){fRejFactorBkgUpgrade=value;}
  void SetIsKeepOnlyCdeuteronSignal(Bool_t isSig){fIsKeepOnlyCdeuteronSignal=isSig;}
  void SetNSoftPionRotations(Int_t nrot){nrot < 0 ? Printf("Cannot set negative number of rotations, setting 0"), fNRotations=0 : fNRotations=nrot;}
  void SetMinAndMaxRotationAngles(Double_t minRot,Double_t maxRot){fMinAngleForRot=minRot;fMaxAngleForRot=maxRot;}
  void SetPDGcodeForFiducialYreco(Int_t pdgcode){fPdgFiducialYreco=pdgcode;}

  // set mass range for ttree filling
  void SetMassRangeTTreeFill( Double_t min, Double_t max )  {flowMass_treeFill=min;  fhighMass_tree_Fill=max;}

  // Sc peak studies in MC
  void SetStudyScPeakMC(Bool_t flag)  {fStudyScPeakMC=flag;}

  // Change the min pT for the soft pion
  void SetMinPtSoftPion(Double_t pTmin) {fMinPtSoftPion=pTmin;}

  void SetNSigmaPreFilterPID(Double_t nsigma) {fNSigmaPreFilterPID=TMath::Abs(nsigma);}

  // switch on/off the ev. sel ev. selection (useful to run on ITS2-ITS3 upgrade MC)
  void SetApplyEvSel(Bool_t flag){fApplyEvSel=flag;}

/*   void SetDoMCAcceptanceHistos(Bool_t doMCAcc=kTRUE){fStepMCAcc=doMCAcc;} */
/*   void SetCutOnDistr(Bool_t cutondistr=kFALSE){fCutOnDistr=cutondistr;} */
/*   void SetUsePid4Distr(Bool_t usepid=kTRUE){fUsePid4Distr=usepid;} */
/*   void SetFillOnlyD0D0bar(Int_t flagfill){fFillOnlyD0D0bar=flagfill;} */
/*   void SetFillVarHists(Bool_t flag) {fFillVarHists=flag;} */
/*   void SetFillPtHistos(Bool_t flag) {fFillPtHist=flag;} */
/*   void SetFillYHistos(Bool_t flag) {fFillYHist=flag;} */
/*   void SetFillImpactParameterHistos(Bool_t flag) {fFillImpParHist=flag;} */
/*   void SetFillSparses(Bool_t flag) {fFillSparses=flag;} */
/*   void SetRejectSDDClusters(Bool_t flag) { fIsRejectSDDClusters=flag; } */
/*   void SetUseSelectionBit(Bool_t flag) { fUseSelectionBit=flag; } */
/*   void SetWriteVariableTree(Bool_t flag) { fWriteVariableTree=flag; } */
/*   void SetDrawDetSignal(Bool_t flag) { fDrawDetSignal=flag; } */
/*   void SetPIDCheck(Bool_t flag) { fPIDCheck=flag; } */
/*   void SetUseQuarkLevelTag(Bool_t opt){fUseQuarkTagInKine=opt;} */
/*   void SetPileupRejectionVZEROTPCout(Bool_t flag) {fEnablePileupRejVZEROTPCout=flag;} */
/*   void SetFillSubSampleHist(Bool_t flag) {fFillSubSampleHist=flag;} */


/*   Bool_t GetCutOnDistr() const {return fCutOnDistr;} */
/*   Bool_t GetUsePid4Distr() const {return fUsePid4Distr;} */
/*   Int_t  GetFillOnlyD0D0bar() const {return fFillOnlyD0D0bar;} */
/*   Bool_t GetFillVarHists() const {return fFillVarHists;} */
/*   Bool_t GetFillPtHistos() const {return fFillPtHist;} */
/*   Bool_t GetFillYHistos() const {return fFillYHist;} */
/*   Bool_t GetFillImpactParameterHistos() const {return fFillImpParHist;} */
/*   Int_t  GetSystem() const {return fSys;} */
/*   Bool_t GetRejectSDDClusters() const { return fIsRejectSDDClusters; } */
/*   Bool_t GetUseSelectionBit() const { return fUseSelectionBit; } */
/*   Bool_t GetWriteVariableTree() const {return fWriteVariableTree;} */
/*   Bool_t GetDrawDetSignal() const {return fDrawDetSignal;} */
/*   Bool_t GetPIDCheck() const {return fPIDCheck;} */
/*   Bool_t GetFillSubSampleHist() const {return fFillSubSampleHist;} */

 private:

  AliAnalysisTaskSEXicTopKpi(const AliAnalysisTaskSEXicTopKpi &source);
  AliAnalysisTaskSEXicTopKpi& operator=(const AliAnalysisTaskSEXicTopKpi& source); 
/*   void	   DrawDetSignal(AliAODRecoDecayHF2Prong *part, TList *ListDetSignal); */

/*   void     FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAODMCHeader *mcHeader, AliRDHFCutsD0toKpi *cuts, TList *listout); */
/*   void     FillVarHists(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout); */
/*   AliAODVertex* GetPrimaryVtxSkipped(AliAODEvent *aodev); */
/*   void CreateImpactParameterHistos(); */
/*   Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const; */
/*   Float_t GetTrueImpactParameter(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partD0) const ; */

  // calculate weight to treat reco true Lc as Xic (mfaggin)
  
  void SigmaCloop(AliAODRecoDecayHF3Prong *io3Prong,AliAODEvent *aod,Int_t massHypothesis,Double_t mass1, Double_t mass2,Double_t *pointS,Int_t resp_onlyPID,Bool_t *arrayPIDselpKpi=0x0,Bool_t *arrayPIDselpiKpi=0x0,Int_t itrack1=-1,Int_t itrack2=-1,Int_t itrackThird=-1,AliAODMCParticle *pSigmaC=0x0,Int_t checkorigin=-1,Int_t decay_channel=0);
  void FillArrayVariableSparse(AliAODRecoDecayHF3Prong *io3Prong,AliAODEvent *aod,Double_t *point,Int_t massHypothesis);  
  Double_t Weight_fromLc_toXic(AliAODMCParticle* p, AliAODMCParticle* prong);
  void PrepareTracks(AliAODEvent *aod,TClonesArray *mcArray=0x0, AliAODMCHeader *mcHeader = 0x0);
  Int_t ConvertXicMCinfo(Int_t infoMC);
  AliAODMCParticle* MatchRecoCandtoMC(AliAODRecoDecayHF3Prong *io3Prong,Int_t &isTrueLambdaCorXic,Int_t &checkOrigin);
  AliAODMCParticle* MatchRecoCandtoMCAcc(AliAODRecoDecayHF3Prong *io3Prong,Int_t &isTrueLambdaCorXic,Int_t &checkOrigin);
  void LoopOverGenParticles();
  void LoopOverFilteredCandidates(TClonesArray *lcArray,AliAODEvent *aod);  
  AliAnalysisVertexingHF *fvHF;   //!<! temporary object for filling reco cands
  AliRDHFCutsD0toKpi *fCuts;      //  Cuts 
  //AliRDHFCutsLctopKpi *fCutsLc;  // Lc Cuts
  AliRDHFCutsXictopKpi *fCutsXic;  // Xic Cuts
  AliNormalizationCounter *fCounter;//!<! AliNormalizationCounter on output slot 5
  AliPIDResponse *fPidResponse; //!<!PID response 
  TClonesArray *fmcArray;        //!<!pointer to mc array
  Bool_t    fReadMC;              ///  flag for MC array: kTRUE = read it, kFALSE = do not read it
  Int_t     fAnalysisType;        /// defines loops and particle targets: 0=default (=all)  1= Xic->pKpi, 2=Lc->pKpi + displacement, 3= Sigma_c, 4= Xicc->Xic pi
  Bool_t fRecalPrimVtx;                  /// switch on/off recalculation of prim vtx: note that in pp and p-Pb it will be set off by default
  TH1F     *fNentries;            //!<! histogram with number of events on output slot 3
  TH1F *fhistMonitoring;         //!<! monitoring histogram
  TArrayI *ftrackArraySel;  //!<! array of selected tracks for internal use
  TArrayI *ftrackArraySelSoftPi;  //!<! array of selected tracks with soft pion cuts, for internal use
  Int_t fnSel;  //!<! number of selected tracks 
  Int_t fnSelSoftPi;  //!<! number of selected tracks with soft-pion cuts
  TArrayI *ftrackSelStatusProton; //!<! array with flags
  TArrayI *ftrackSelStatusKaon;//!<! array with flags
  TArrayI *ftrackSelStatusPion;//!<! array with flags
  Int_t     fSys;                 /// fSys=0 -> p-p; fSys=1 ->pPb; fSys=2 ->PbPb (PV vtx is recalculated only in pp)
  Int_t     fAODProtection;       /// flag to activate protection against AOD-dAOD mismatch.
                                  /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  Int_t fLikeSign;                /// flag for like sign
  AliESDtrackCuts *fESDtrackCutsProton; //
  AliESDtrackCuts *fESDtrackCutsKaon;//
  AliESDtrackCuts *fESDtrackCutsPion;//
  AliESDtrackCuts *fESDtrackCutsSoftPion;//
  AliAODVertex *fprimVtx;//! pointer to prim. vertex
  TH2F *fhistInvMassCheck;//! hist with generic inv. mass distr (for checks)

  //TH3F *fhistMCSpectrumAccLc;//! hist with MC spectrum of cand in acceptance
  THnF *fhistMCSpectrumAccLc;//! hist with MC spectrum of cand in acceptance

  THnSparseF *fhistMCSpectrumAccLcFromSc;//! hist with MC spectrum of cand in acceptance

  //TH3F *fhistMCSpectrumAccSc;//! hist with MC spectrum of cand in acceptance
  THnF *fhistMCSpectrumAccSc;//! hist with MC spectrum of cand in acceptance

  //TH3F *fhistMCSpectrumAccXic;//! hist with MC spectrum of cand in acceptance
  THnF *fhistMCSpectrumAccXic;//! hist with MC spectrum of cand in acceptance

  TH2F *fhistMCSpectrumAccCdeuteron;//! hist with MC spectrum of cand in acceptance
  THnSparseF* fhSparseAnalysis;//! sparse for analysis
  THnSparseF* fhSparseAnalysisSigma;//! sparse for analysis of SigmaC (with deltaM)
  THnSparseF* fhSparsePartReco;//! sparse for single track efficiency (reco spectra)
  THnSparseF* fhSparsePartGen;//! sparse for single track efficiency (gen spectra)
  TH1F* fCosPointDistrAll; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fCosPointDistrAllFilter; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fCosPointDistrSignal; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fCosPointDistrSignalFilter; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fDist12Signal; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fDist12SignalFilter; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fDist12All; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fDist12AllFilter; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fDist23Signal; //!<! histo storing variable for debugging (pt integrted)
  TH1F* fDist23All; //!<! histo storing variable for debugging (pt integrted)
  TH1F *fDist23AllFilter;//!<! histo storing variable for debugging (pt integrted)
  TH2F *fVtxResXPt; //!<! histo for vertex resolution in X vs pT
  TH2F *fVtxResYPt; //!<! histo for vertex resolution in Y vs pT
  TH2F *fVtxResZPt; //!<! histo for vertex resolution in Z vs pT
  TH2F *fVtxResXYPt; //!<! histo for vertex resolution in XY vs pT
  TH2F *fVtxResXYZPt; //!<! histo for vertex resolution in XYZ vs pT
  TH2F *fPrimVtxResXPt; //!<! histo for primary vertex resolution in X vs pT
  TH2F *fPrimVtxResYPt; //!<! histo for primary vertex resolution in Y vs pT
  TH2F *fPrimVtxResZPt; //!<! histo for primary vertex resolution in Z vs pT
  TH2F *fDecayLResXPt; //!<! histo for decay length resolution in X vs pT
  TH2F *fDecayLResYPt; //!<! histo for decay length resolution in Y vs pT
  TH2F *fDecayLResZPt; //!<! histo for decay length resolution in Z vs pT
  TH2F *fDecayLResXYPt; //!<! histo for decay length resolution in XY vs pT
  TH2F *fDecayLResXYZPt; //!<! histo for decay length resolution in XYZ vs pT
  TH2F *fnSigmaPIDtofProton; //!<! histo for monitoring PID performance
  TH2F *fnSigmaPIDtofPion; //!<! histo for monitoring PID performance
  TH2F *fnSigmaPIDtofKaon; //!<! histo for monitoring PID performance
  TH2F *fnSigmaPIDtpcProton; //!<! histo for monitoring PID performance
  TH2F *fnSigmaPIDtpcPion; //!<! histo for monitoring PID performance
  TH2F *fnSigmaPIDtpcKaon; //!<! histo for monitoring PID performance
  TH2F *fProtonID; //!<! histo for purity of PID
  TH2F *fKaonID; //!<! histo for purity of PID
  TH2F *fPionID; //!<! histo for purity of PID
  TList *fOutput;//! Output List
  AliVertexerTracks *fVertexerTracks;//!<! vertexer
  Bool_t fSetTrackCutLcFilteringPP; /// flag to force esd track cuts used for Lc filtering
  Int_t fCutSelLevel; /// flag to define cuts used online
  Bool_t fApplykFirst;/// flag to apply kFirst selection at track level for pt<fMaxPtTrackkFirst (needed just to avoid pt calculations if not needed)
  Double_t fMaxPtTrackkFirst;//
  Double_t fMaxVtxChi2Cut;//
  //Int_t fFillTree;//
  Bool_t fFillTree;//
  TTree *fTreeVar;//

/*   TList    *fOutputMass;          //!<! list send on output slot 1 */
/*   TList    *fOutputMassPt;        //!<! list send on output slot 6 */
/*   TList    *fOutputMassY;        //!<! list send on output slot 9 */
/*   TList    *fDistr;               //!<! list send on output slot 2 */

/*   THnSparseF *fMCAccPrompt;       //!<!histo for StepMCAcc for D0 prompt (pt,y,ptB) */
/*   THnSparseF *fMCAccBFeed;        //!<!histo for StepMCAcc for D0 FD (pt,y,ptB) */
/*   Bool_t fStepMCAcc;              // flag to activate histos for StepMCAcc */

/*   THnSparseF *fHistMassPtImpParTC[5];   //!<! histograms for impact paramter studies */
/*   Int_t     fArray;               ///  can be D0 or Like Sign candidates */

/*   Bool_t    fCutOnDistr;          ///  flag to decide if apply cut also on distributions: 0 no cuts, 1 looser cuts, 2 tighter/ cuts */
/*   Bool_t    fUsePid4Distr;        //  flag to use the particle identification to fill the signal histograms of distributions. It has effect only with fReadMC=kFALSE */

/*   Int_t     fNPtBins;             ///  number of pt bins */
/*   Double_t  fLsNormalization;     ///  normalization */
/*   Int_t     fFillOnlyD0D0bar;     /// flag to fill mass histogram with D0/D0bar only (0 = fill with both, 1 = fill with D0 only, 2 = fill with D0bar only) */
/*   TObjArray fDaughterTracks;      /// keeps the daughter tracks */
/*   Int_t     fIsSelectedCandidate; /// selection outcome */
/*   Bool_t    fFillVarHists;        /// flag to enable filling variable histos */

/*   Bool_t    fIsRejectSDDClusters; /// flag to reject events with SDD clusters */
/*   Bool_t    fFillPtHist;          /// flag to fill Pt and Impact Parameter Histograms */
/*   Bool_t    fFillYHist;          /// flag to fill Y Histograms */
/*   Bool_t    fFillImpParHist;      /// flag to fill Pt and Impact Parameter Histograms */
/*   Bool_t    fFillSubSampleHist;    /// flag to fill SubSample histogram */
/*   Int_t     fEventCounter; /// event counter used for sub sample test */
/*   Bool_t    fUseSelectionBit;     /// flag to check or not the selection bit */

/*   Bool_t    fWriteVariableTree;       /// flag to decide whether to write the candidate variables on a tree variables */
/*   TTree    *fVariablesTree;           //!<! tree of the candidate variables after track selection on output slot 7 */
/*   Double_t *fCandidateVariables;      //!<!  variables to be written to the tree */
/*   Bool_t	fPIDCheck;			/// flag to decide whether to fill "PID = x" bins in fNentrie */
/*   Bool_t    fDrawDetSignal;		/// flag to decide whether to draw the TPC dE/dx and TOF signal before/after PID */
/*   Bool_t fUseQuarkTagInKine;            // flag for quark/hadron level identification of prompt and feeddown */
/*   Bool_t fFillSparses;                  // flag to activate THnSparse  */
/*   THnSparseF *fhStudyImpParSingleTrackSign; //!<! sparse with imp par residual cuts for MC */
/*   THnSparseF *fhStudyImpParSingleTrackCand;  //!<! sparse with imp par residual cuts for Data */
/*   THnSparseF *fhStudyImpParSingleTrackFd;   //!<! sparse with imp par residual cuts for MC */
/*   TList	   *fDetSignal;		//!<!Detector signal histograms (on output slot 8) */
/*   TH2F *fhMultVZEROTPCoutTrackCorrNoCut;  //!<! */
/*   TH2F *fhMultVZEROTPCoutTrackCorr;  //!<! */
/*   Bool_t    fEnablePileupRejVZEROTPCout; */

  // downsampling for fTreeVar filling
  Float_t fpT_down;         /// pT value that discriminates between low and high pT for the downsampling 
  Float_t fLowpT_down;      /// downsampling factor at low pT
  Float_t fHighpT_down;     /// downsampling factor at high pT
  // pT limits for TTree filling
   Float_t fminpT_treeFill;   /// min. pT
   Float_t fmaxpT_treeFill;   /// max. pT
  Bool_t fCompute_dist12_dist23;  /// flag to require the calculation of dist12 and dist23

  Bool_t fExplore_PIDstdCuts;     /// flag to switch on the exporation of PID cuts with standard strategy
  Bool_t fOnlyBayesPIDbin;        /// dirty solution: flag to reduce the axes in the reco sparses ---> make the merging easier (mfaggin)
  Bool_t fExplPID_BayesOnlyProt;  /// flag to include the PID selection with Bayes approach only for proton

  Double_t fLcMassWindowForSigmaC; /// lc mass window for used in sigma_C loop
  Double_t fSigmaCDeltaMassWindow; /// mass window for accetping sigma_C candidate
  Bool_t fSigmaCfromLcOnTheFly; /// switch to use on-the-fly Lc or filtered Lc from delta file
  Bool_t fCheckOnlyTrackEfficiency;// flag for filling only the single-track sparse and return
  Bool_t fIsCdeuteronAnalysis;// flag for doing the c deuteron analysis (inv mass)
  Bool_t fIsKeepOnlyCdeuteronSignal;// flag for keeping only c deuteron signal
  Bool_t fIsXicUpgradeAnalysis; // flag for ITS2-3 studies for Xic
  Bool_t fIsKeepOnlySigXicUpgradeAnalysis; // flag for keeping only Xic signal
  Bool_t fIsKeepOnlyBkgXicUpgradeAnalysis; // flag for keeping only Xic bkg
  Double_t fRejFactorBkgUpgrade;  // rejection factor for background reconstruction in upgrade studies
  Int_t fNRotations;    // number of rotations performed on soft pion, to study SigmaC background shape; 0 = no rotations, 1 -> single rotations by fMinAngleForRot, 2 -> fNRotations from fMinAngleForRot to fMaxAngleForRot
  Double_t fMinAngleForRot;//
  Double_t fMaxAngleForRot;//
  Int_t   fPdgFiducialYreco;// pdg code of particle that will be used to compute fiducial acceptance at reco level

  // mass ranges to store candidates in the TTree
  Double_t flowMass_treeFill;
  Double_t fhighMass_tree_Fill;

  // THnSparse for Sc peak studies in MC
  Bool_t fStudyScPeakMC;
  THnSparseF* fhsparseMC_ScPeak; //!

  // double to change the min pT for the soft pion
  Double_t fMinPtSoftPion;  // !

  // histogram to look at the vtx_z in MC
  TH1D* fZvtx_gen_within10_MC; //!<!
  TH1D* fZvtx_gen_noSel10_MC; //!<!
  TH1D* fZvtx_reco_noSel10_MC; //!<!

  // histogram to count candidates
  TH1F* fCandCounter; //!<!
  TH1F* fCandCounter_onTheFly; //!<!

  // pT distribution of candidate soft pion tracks
  TH1F* fPtSoftPionCand; //!<!
  // pT distribution of candidate soft pion tracks inside SigmaC loop
  TH1F* fPtSoftPionCand_insideScLoop; //!<!

  Double_t fNSigmaPreFilterPID; // number of sigma for TPC and TOF pre-filtering PID on tracks
  
  // bool to remove ev. selection (useful to run on ITS2-ITS3 upgrade MC)
  Bool_t fApplyEvSel;

  // bool to keep only the Bayes PID- based PID selections
  Bool_t fNoStdPIDcases;

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEXicTopKpi,18); /// AliAnalysisTaskSE for Xic->pKpi
  /// \endcond
};

#endif

