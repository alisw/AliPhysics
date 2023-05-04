#ifndef AliAnalysisTaskSELcTopK0sCorrelations_H
#define AliAnalysisTaskSELcTopK0sCorrelations_H

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskSELcTopK0sCorrelations.h $ */

//*************************************************************************
// Class AliAnalysisTaskSELcTopK0sCorrelations
// AliAnalysisTaskSE for Lambdac candidates (2Prongs) and hadrons correlations
// Authors:
// A.Palasciano,  antonio.palasciano@ba.infn.it 
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <TProfile.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliHFMLResponseLambdactopK0s.h"
#include "AliHFMLVarHandler.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliHFCorrelator.h"
#include "AliNormalizationCounter.h"
#include "AliHFOfflineCorrelator.h"
#include "AliD0hCutOptim.h"
#include "AliAnalysisTaskSELambdac.h"
#include "AliAODRecoCascadeHF.h"

using std::vector;

class AliAODEvent;

class AliAnalysisTaskSELcTopK0sCorrelations : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSELcTopK0sCorrelations();
  AliAnalysisTaskSELcTopK0sCorrelations(const char *name,AliRDHFCutsLctoV0* cutsLambdac, AliHFAssociatedTrackCuts* cutsTrk);
  virtual ~AliAnalysisTaskSELcTopK0sCorrelations();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  enum PartType {kTrack,kKCharg,kK0};
  enum FillType {kSE, kME}; //for single event or event mixing histos fill
  //enum TreeFill {kNoTrees, kFillTrees, kFillCutOptTree};
  //enum SpeedType {kAllBins, kOneBinSB, kOneBinSBandS};
  
  enum {kNoTrees, kFillTrees, kFillCutOptTree};
  enum {kAllBins, kOneBinSB, kOneBinSBandS};

  void SetReadMC(Bool_t readMC=kFALSE){fReadMC=readMC;}
  void SetMCReconstructedTracks(Bool_t recoTrMC=kTRUE){fRecoTr=recoTrMC;}
  void SetMCReconstructedLambdac(Bool_t recoLambdacMC=kTRUE){fRecoLambdac=recoLambdacMC;}
  void SetMCSelEventType(Bool_t sel=kFALSE){fSelEvType=sel;}
  void SetFillOnlyLambdacLambdacbar(Int_t flagfill){fFillOnlyLambdacLambdacbar=flagfill;}
  void SetSystem(Int_t sys){fSys=sys;}
  void SetRejectSDDClusters(Bool_t flag) {fIsRejectSDDClusters=flag; }
  void SetFillGlobalPlots(Bool_t fill=kTRUE){fFillGlobal=fill;}
  void SetSoftPiFlag(Bool_t piflag) {fSoftPiCut=piflag;}
  void SetMEAxisThresh(Bool_t methresh) {fMEAxisThresh=methresh;}
  void SetKaonCorrelations(Bool_t kaonCorr) {fKaonCorr=kaonCorr;}
  void SetAODMismatchProtection(Int_t opt=0) {fAODProtection=opt;}
  void SetPurityStudies(Bool_t puritystudies=kFALSE) {fPurityStudies=puritystudies;}

  Int_t  GetReadMC() const {return fReadMC;}
  Int_t  GetMCReconstructedTracks() const {return fRecoTr;}
  Int_t  GetMCReconstructedLambdac() const {return fRecoLambdac;}
  Int_t  GetMCSelEventType() const {return fSelEvType;}
  Int_t  GetFillOnlyLambdacLambdacbar() const {return fFillOnlyLambdacLambdacbar;}
  Int_t  GetSystem() const {return fSys;}
  Bool_t GetRejectSDDClusters() const {return fIsRejectSDDClusters;}
  Bool_t GetFillGlobalPlots() const {return fFillGlobal;}
  Double_t GetEtaForCorrel() {return fEtaForCorrel;}
  Double_t GetMultEv() {return fMultEv;}
  Bool_t GetSoftPiFlag() const {return fSoftPiCut;}
  Bool_t GetMEAxisThresh() const {return fMEAxisThresh;}
  Bool_t GetKaonCorrelations() const {return fKaonCorr;}
  Bool_t GetFillTrees() const {return fFillTrees;}

  //correlations setters/printers
  void SetNPtBinsCorr(Int_t nbins) {fNPtBinsCorr = nbins;}
  void SetPtBinsLimsCorr(Double_t* ptlims) {for(int i=0;i<=fNPtBinsCorr;i++) {fBinLimsCorr.push_back(ptlims[i]);}}
  void SetPtBinsLimsCorr(Float_t* ptlims) {for(int i=0;i<=fNPtBinsCorr;i++) {fBinLimsCorr.push_back((Double_t)ptlims[i]);}}
  void SetPtTreshLow(Double_t* pttreshlow) {for(int i=0;i<fNPtBinsCorr;i++) {fPtThreshLow.push_back(pttreshlow[i]);}}
  void SetPtTreshUp(Double_t* pttreshup) {for(int i=0;i<fNPtBinsCorr;i++) {fPtThreshUp.push_back(pttreshup[i]);}}
  void SetLSBLowLim(Double_t* LSBLowLim) {for(int i=0;i<fNPtBinsCorr;i++) {fLSBLowLim.push_back(LSBLowLim[i]);}}
  void SetLSBHighLim(Double_t* LSBUppLim) {for(int i=0;i<fNPtBinsCorr;i++) {fLSBUppLim.push_back(LSBUppLim[i]);}}  
  void SetRSBLowLim(Double_t* RSBLowLim) {for(int i=0;i<fNPtBinsCorr;i++) {fRSBLowLim.push_back(RSBLowLim[i]);}}
  void SetRSBHighLim(Double_t* RSBUppLim) {for(int i=0;i<fNPtBinsCorr;i++) {fRSBUppLim.push_back(RSBUppLim[i]);}}
  void SetSignLowLim(Double_t* SignLowLim) {for(int i=0;i<fNPtBinsCorr;i++) {fSignLowLim.push_back(SignLowLim[i]);}}
  void SetSignHighLim(Double_t* SignUppLim) {for(int i=0;i<fNPtBinsCorr;i++) {fSignUppLim.push_back(SignUppLim[i]);}}  
  void SetLeftSignReg_LowPt(Double_t leftlow) {fSignLeft_LowPt=leftlow;}
  void SetRightSignReg_LowPt(Double_t rightlow) {fSignRight_LowPt=rightlow;}
  void SetLeftSignReg_HighPt(Double_t lefthigh) {fSignLeft_HighPt=lefthigh;}
  void SetRightSignReg_HighPt(Double_t righthigh) {fSignRight_HighPt=righthigh;}
  void SetPtAssocLim(Double_t pTlim) {fPtAssocLimit=pTlim;}
  
  void PrintBinsAndLimits();
  Int_t PtBinCorr(Double_t pt) const;
  void SetEvMixing(Bool_t mix) {fMixing=mix;}
  void SetEtaForCorrel(Double_t etacorr) {fEtaForCorrel=etacorr;}
  void SetSpeed(Int_t speed) {fSpeed=speed;}
  void SetMergePools(Bool_t mergepools) {fMergePools=mergepools;}
  void SetUseLceff(Bool_t UseLceff) {fUseLceff=UseLceff;}
  void SetUseTrackeff(Bool_t useTrackeff) {fUseTrackeff=useTrackeff;}
  void SetMinDPt(Double_t minDPt) {fMinDPt=minDPt;}
  void SetFillTrees(Int_t fillTrees, Double_t fractAccME) {fFillTrees=fillTrees; fFractAccME=fractAccME;}
  void SetCentralityV0(Double_t V0min, Double_t V0max) {fV0CentMin=V0min; fV0CentMax=V0max;}  
  void SetTrackletRange(Double_t trkmin, Double_t trkmax) {fTrkMultMin=trkmin; fTrkMultMax=trkmax;}
  void SetAnalysisVsMult(Bool_t v2anal) {fVsMultAnalysis=v2anal;}
 
  void SetUseNtrklWeight(Bool_t flag=kTRUE) {fUseNtrklWeight=flag;}
  void SetHistNtrklWeight(TH1D* h) {
    if(fHistNtrklWeight) delete fHistNtrklWeight;
    fHistNtrklWeight = new TH1D(*h);
  }

  void SetEqualizeTracklets(Bool_t flag) {fEqualizeTracklets=flag;}  
  void SetReferenceMultiplicity(Double_t refmult) {fRefMult=refmult;}
  void SetMultiplVsZProfile(TProfile* hprof, Int_t index){
    if(fTrackletProfiles[index]) delete fTrackletProfiles[index];
    fTrackletProfiles[index]=new TProfile(*hprof);
  }
  
  /// methods for ML application
  void SetDoMLApplication(Bool_t flag = kTRUE, Bool_t isMultiClass = kFALSE) {fApplyML = flag; fMultiClassML = isMultiClass;}
  void SetMLConfigFile(TString path = "") {fConfigPath = path;}
  void SetMLBinsForSparse(Int_t nbins = 300, Double_t min = 0.85, Double_t max = 1.){ fNMLBins[0] = nbins; fMLOutputMin[0] = min; fMLOutputMax[0] = max;}
  void SetMultiClassMLBinsForSparse(Int_t nbinsBkg = 100,
                                    Int_t nbinsPrompt = 100,
                                    Int_t nbinsFD = 100,
                                    Double_t minBkg = 0., Double_t maxBkg = 1.,
                                    Double_t minPrompt = 0., Double_t maxPrompt = 1.,
                                    Double_t minFD = 0., Double_t maxFD = 1.) 
  {
    fNMLBins[0] = nbinsBkg; fNMLBins[1] = nbinsPrompt; fNMLBins[2] = nbinsFD;
    fMLOutputMin[0] = minBkg; fMLOutputMin[1] = minPrompt; fMLOutputMin[2] = minFD;
    fMLOutputMax[0] = maxBkg; fMLOutputMax[1] = maxPrompt; fMLOutputMax[2] = maxFD;
  }
  void SetMinimalVarForMLSparse(Bool_t flag = kTRUE) {fUseMinimalVarForSparse = flag;}
  /// methods for ML tree creation
  void SetCreateMLTree(Bool_t flag = kTRUE) {fCreateMLtree = flag;}
  void SetMLTreePIDopt(int opt) {fPIDopt = opt;} // default AliHFMLVarHandlerLctopKpi::kNsigmaDetAndCombPID
  void SetMLTreeAddTrackVar(Bool_t flag = kTRUE) {fAddSingleTrackVar = flag;}
  void SetMLTreeAddImpParProd(Bool_t flag = kTRUE) {fAddImpParProdProngs = flag;}
  void SetFillOnlySignalInMLtree(Bool_t opt = kTRUE) {
    if(fReadMC) fFillOnlySignal = opt;
    else {
      if(opt)
        AliError("fReadMC has to be kTRUE");
    }
  }
  void SetMLTreeAddNtracklets(Bool_t flag = kTRUE) {fAddNtracklets = flag;}

  void EnableMLTreeEvtSampling(Float_t fractokeep, ULong_t seed) {
    fEnableEvtSampling = kTRUE;
    fFracEvtToKeep = fractokeep;
    fSeedSampling = seed;
  }
  void EnableMLTreeCandSampling(Float_t fractokeep, Float_t maxptsampling) {
    fEnableCandSampling = kTRUE;
    fFracCandToKeep = fractokeep;
    fMaxCandPtSampling = maxptsampling;
  }


 private:

  AliAnalysisTaskSELcTopK0sCorrelations(const AliAnalysisTaskSELcTopK0sCorrelations &source);
  AliAnalysisTaskSELcTopK0sCorrelations& operator=(const AliAnalysisTaskSELcTopK0sCorrelations& source); 
  void FillMassHists(AliAODRecoCascadeHF *part, TClonesArray *arrMC, AliRDHFCutsLctoV0 *cuts, TList *listout, AliAODEvent *aod);
  Int_t CheckLambdacOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  //correlation methods
  void CreateCorrelationsObjs();
  void CalculateCorrelations(AliAODRecoCascadeHF* d, Int_t labLambdac=-1, TClonesArray* mcArray=0x0);
  void CalculateCorrelationsMCKine(AliAODMCParticle* d, TClonesArray* mcArray=0x0);
  void FillSparsePlots(TClonesArray* arrayMC, Double_t mInv[], Int_t origLambdac, Int_t PdgLambdac, AliReducedParticle* track, Int_t ptbin, Int_t type, Int_t softpiME, Double_t wg=1.);
  Int_t CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  Bool_t IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* mcArray) const;
  void FillTreeLambdac(AliAODRecoCascadeHF* d, AliAODEvent* aod);  
  void FillTreeTracks(AliAODEvent* aod);  
  void FillTreeLambdacForCutOptim(AliAODRecoCascadeHF* d, AliAODEvent* aod);  
  void ResetBranchD();
  void ResetBranchTracks();
  void ResetBranchDForCutOptim();
  Bool_t AcceptTrackForMEOffline(Double_t pt);
  void FillPurityPlots(TClonesArray* mcArray, AliReducedParticle* track, Int_t ptbin, Double_t deltaphi);
  Double_t GetNtrklWeight(Int_t ntrkl); 
  TProfile* GetEstimatorHistogram(const AliVEvent* event);
  bool ReconstructKFLc(AliAODRecoCascadeHF *cand);
  bool CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau);
  
  Int_t             	   fNPtBinsCorr;        // number of pt bins per correlations
  std::vector<Double_t>  fBinLimsCorr;        // limits of pt bins per correlations
  std::vector<Double_t>  fPtThreshLow;        // pT threshold of hadrons - low
  std::vector<Double_t>  fPtThreshUp;         // pT threshold of hadrons - up
  std::vector<Double_t>  fLSBLowLim;          // Left SB lower lim
  std::vector<Double_t>  fLSBUppLim;          // Left SB lower lim
  std::vector<Double_t>  fRSBLowLim;          // Right SB upper lim
  std::vector<Double_t>  fRSBUppLim;          // Right SB upper lim
  std::vector<Double_t>  fSignLowLim;         // Left signal region lower lim (for fSpeed==2)
  std::vector<Double_t>  fSignUppLim;         // Left signal region lower lim (for fSpeed==2)
  std::vector<Int_t>     fDaughTrackID;       // ID of tagged daughters
  std::vector<Int_t>     fDaughTrigNum;	      // ID of D-trigger for daughters	

  Int_t     fEvents;		  	// EventCounter
  Bool_t    fAlreadyFilled;	  	// Lambdac in an event already analyzed (for track distribution plots)
  Int_t	    fNtrigD;			// counter on number of D triggers filled (for association with daughter tracks in TTrees)
  TList    *fOutputMass;          	//!list send on output slot 1
  TList    *fOutputCorr;	  	//!list of correlation histos, output slot 5
  TList    *fOutputStudy;	  	//!list of histos with MC distributions, output slot 6
  TH1F     *fNentries;            	//!histogram with number of events on output slot 2
  AliRDHFCutsLctoV0 *fCutsLambdac;    	// Cuts for Lambdac, output 3
  AliHFAssociatedTrackCuts *fCutsTracks;// Cuts for tracks and K0, output 7
  AliHFCorrelator* fCorrelatorTr;	// Correlator for tracks
  AliHFCorrelator* fCorrelatorKc;	// Correlator for charged K
  AliHFCorrelator* fCorrelatorK0;	// Correlator for K0
  Bool_t    fReadMC;              	// flag for MC array: kTRUE = read it, kFALSE = do not read it
  Bool_t    fRecoTr;   		       	// flag for using MC reconstructed (kTRUE) or pure kinematic MC (kFALSE) - Associated tracks
  Bool_t    fRecoLambdac;   		       	// flag for using MC reconstructed (kTRUE) or pure kinematic MC (kFALSE) - Lambdac
  Bool_t    fSelEvType;		       	// flag for enabling selection of event tpye (PP, GS, FE, ...) on MC analysis
  Bool_t    fMixing;			// flag to enable also event mixing
  AliNormalizationCounter *fCounter;	//!AliNormalizationCounter on output slot 4
  Int_t     fNPtBins;             	// Number of pt bins
  Int_t     fFillOnlyLambdacLambdacbar;     	// flag to fill mass histogram with Lambdac/Lambdacbar only (0 = fill with both, 1 = fill with Lambdac only, 2 = fill with Lambdacbar only)
  Int_t     fIsSelectedCandidate; 	// selection outcome
  Int_t     fSys;                 	// fSys=0 -> p-p; fSys=1 ->PbPb
  Double_t  fEtaForCorrel;		// cut for Lambdac eta to enable correlation with associated particles
  Bool_t    fIsRejectSDDClusters; 	// flag to reject events with SDD clusters
  Bool_t    fFillGlobal;          	// flag to fill global plots (in loops on tracks and V0 for each event)
  Double_t  fMultEv;			// event multiplicity (for trigger eff), if in terms of tracklets, is equalized (if asked!)!
  Double_t  fMultEvOrig;      // event multiplicity (for trigger eff), not equalized!
  Double_t  fMultEvV0M;     // event multiplicity (for trigger eff)
  Double_t  fMultEvV0MEqual;     // event multiplicity (for trigger eff)
  Double_t  fCentEvV0M;     // event multiplicity (for trigger eff)
  Double_t  fzVtx;				// event zVtx position (for track eff)
  Bool_t    fSoftPiCut;			// flag to activate soft pion cut on Data
  Bool_t    fMEAxisThresh;		// flag to fill threshold axis in ME plots
  Bool_t    fKaonCorr;			// enables correlations of Lambdac-Kcharg and Lambdac-K0
  Double_t  fSignLeft_LowPt;		// Left bound of "signal region" range - up to 8 GeV/c
  Double_t  fSignRight_LowPt;		// Right bound of "signal region" range - up to 8 GeV/c
  Double_t  fSignLeft_HighPt;		// Left bound of "signal region" range - from 8 GeV/c
  Double_t  fSignRight_HighPt;		// Right bound of "signal region" range - from 8 GeV/c
  Int_t     fPoolNum;			// Number of the pool for the analyzed event
  //SpeedType fSpeed;			// Speed up the execution removing bins and histos - 0=std, 1=single-SB bins, 2=single-SB and single-S bins
    Int_t fSpeed;			// Speed up the execution removing bins and histos - 0=std, 1=single-SB bins, 2=single-SB and single-S bins
  Bool_t    fMergePools;		// Put all entries from various pools in _pool0 THnSparses (as old approach) - for testing & low stat!
  Bool_t    fUseLceff;			// Use D meson efficiency as weight
  Bool_t    fUseTrackeff;   		// Use track efficiency as weight
  Double_t  fPtAssocLimit;   		// Maximum value for associated pT
  Double_t  fMinDPt;			// Minimum pT of the Lambdac to allow selection
  Double_t  fV0CentMin;         // Minimum V0 centrality for internal event selection (not made by cut object)
  Double_t  fV0CentMax;         // Maximum V0 centrality for internal event selection (not made by cut object)
  Double_t  fTrkMultMin;        // Minimum SPD tracklets in |eta|<1 for internal event selection (not made by cut object)
  Double_t  fTrkMultMax;        // Minimum SPD tracklets in |eta|<1 for internal event selection (not made by cut object)
  Bool_t    fVsMultAnalysis;        // Running v2 in pp analysis

  //TreeFill  fFillTrees;			// Flag to fill ME offline trees
    Int_t  fFillTrees;			// Flag to fill ME offline trees
  Double_t  fFractAccME;		// Fraction of tracks to be accepted in the ME offline
  Int_t     fAODProtection;  	        // flag to activate protection against AOD-dAOD mismatch.
  Bool_t    fPurityStudies;		// flag to activate purity studies (primaries, secondaries, charm and beauth tracks rejected by DCA cut, vs pT and deltaPhi)

  Bool_t    fUseNtrklWeight;            // flag to activate events weighting via Ntracklet distribution data/MC ratio
  TH1D      *fHistNtrklWeight;          // histo with Ntracklets weights
  Double_t  fWeight;                    // Ntrkl weight to apply to events when filling the MC THnSparse (for closure test) and MC purity plots

  Bool_t     fEqualizeTracklets;     //activates the tracklet correction using the TProfiles (data)
  Double_t   fRefMult;               //refrence multiplcity (max of maxes of profiles in dataset)
  TProfile*  fTrackletProfiles[33];  //TProfile with mult vs. Z per period

  AliHFCorrelationBranchD   *fBranchD;
  AliHFCorrelationBranchTr  *fBranchTr;
  AliD0hCutOptim	    *fBranchDCutVars; //for cut optimization!

  TTree	    *fTreeD;			// TTree for ME offline - Lambdac mesons
  TTree	    *fTreeTr;			// TTree for ME offline - Assoc tracks
  TObjArray *fTrackArray;		// Array with selected tracks for association
  Bool_t    fTrackArrayFilled;		// Flag to fill fTrackArray or not (if already filled)
  
  /// variables for ML application
  Bool_t fApplyML = kFALSE;                                           /// flag to enable ML application
  Bool_t fMultiClassML = kFALSE;                                      /// flag to enable multi-class models (bkg, prompt, FD)
  TString fConfigPath = "";                                           /// path to ML config file
  AliHFMLResponseLambdactopK0s* fMLResponse = nullptr;                 //!<! object to handle ML response
  Int_t fNMLBins[3] = {100, 100, 100};                                /// number of bins for ML output axis in THnSparse
  Double_t fMLOutputMin[3] = {0.0, 0.0, 0.0};                         /// min for ML output axis in THnSparse
  Double_t fMLOutputMax[3] = {1.0, 1.0, 1.0};                         /// max for ML output axis in THnSparse

  /// variables for tree creation
  Bool_t fCreateMLtree = kFALSE;
  AliHFMLVarHandler* fMLhandler = nullptr;                //!<! object to handle ML tree creation and filling
  TTree* fMLtree = nullptr;                                           //!<! tree with candidates for ML
  int fPIDopt = AliHFMLVarHandler::kNsigmaDetAndCombPID;  /// option for PID variables
  Bool_t fAddSingleTrackVar = kFALSE;                                 /// option to store single track variables
  Bool_t fAddImpParProdProngs = kFALSE;                               /// option to store LambdacK*Lambdacpi1 and LambdacK*Lambdacpi2 variables
  Bool_t fAddNtracklets = kFALSE;                                     /// option to add Ntracklets 
  Bool_t fFillOnlySignal = kFALSE;                                    /// option to store only signal when using MC
  Bool_t fUseMinimalVarForSparse = kFALSE;                            /// flag to keep only inv. mass, pt and prob. in the sparse

  Bool_t fEnableEvtSampling = kFALSE;                                 /// flag to apply event sampling
  Float_t fFracEvtToKeep = 1.1;                                       /// fraction of events to be kept by event sampling
  ULong_t fSeedSampling = 0;                                          /// seed for sampling
  Bool_t fEnableCandSampling = kFALSE;                                /// flag to apply candidate sampling
  Float_t fFracCandToKeep = 1.1;                                      /// fraction of candidates to be kept by sampling
  Float_t fMaxCandPtSampling = 0.;                                    /// maximun candidate pt to apply sampling

  //KFVariables Lc candidate
  Double_t fKFLcInvMass = 0.;
  Double_t fKFLcPt = 0.;
  Double_t fKFLcEta = 0.;  
  Double_t fKFLcPhi = 0.;  
  Double_t fKFLcRapidity = 0.;

  ClassDef(AliAnalysisTaskSELcTopK0sCorrelations,19); // AliAnalysisTaskSE for Lambdac->Kpi - h correlations
};

#endif

