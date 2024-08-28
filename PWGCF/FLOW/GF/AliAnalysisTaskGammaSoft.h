#ifndef ALIANALYSISTASK_GAMMASOFT__H
#define ALIANALYSISTASK_GAMMASOFT__H
#include "AliAnalysisTaskSE.h"
#include "AliGFW.h"
#include "AliEventCuts.h"
#include "TString.h"
#include "TRandom.h"
#include "TAxis.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "AliPtContainer.h"
#include "AliProfileBS.h"
#include "AliGFWWeights.h"
#include "AliGFWCuts.h"
#include "AliGFWFlowContainer.h"
#include "TF1.h"

class AliAODEvent;
class AliAODTrack;
class AliVEvent;

class AliAnalysisTaskGammaSoft : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskGammaSoft();
  AliAnalysisTaskGammaSoft(const char *name, Bool_t IsMC=kTRUE, TString ContainerSubfix="");
  virtual ~AliAnalysisTaskGammaSoft();
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  Bool_t CheckTrigger(Double_t);
  Bool_t AcceptAOD(AliAODEvent*, Double_t lvtxXYZ[3]);
  Bool_t IsPileupEvent(AliAODEvent* ev, double centrality);
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
  void SetEventCutFlag(Int_t newval) { fEventCutFlag = newval; };
  AliGFW::CorrConfig GetConf(TString head, TString desc, Bool_t ptdif) { return fGFW->GetCorrelatorConfig(desc.Data(),head.Data(),ptdif);};
  void CreateCorrConfigs();
  void LoadCorrectionsFromLists();
  void FillWPCounter(vector<vector<double>> &inarr, double w, double p);
  void SetPtBins(Int_t nBins, Double_t *ptbins);
  void SetEtaBins(Int_t nBins, Double_t *etabins);
  void SetMultiBins(Int_t nBins, Double_t *multibins);
  void SetV0MBins(Int_t nBins, Double_t *multibins);
  void SetDCABins(Int_t nDCABins, Double_t *DCAbins);
  void SetDCABins(Int_t nDCABins, Double_t low, Double_t high);
  void SetEtaMpt(Double_t newval) { fEtaMpt = newval; };
  void SetEtaAcceptance(Double_t newval) { fEtaAcceptance = newval; };
  void SetEtaV2Sep(Double_t newval) { fEtaV2Sep = newval; };
  void SetUseNch(Bool_t newval) { fUseNch = newval; };
  void SetUseUnityParticleWeights(Bool_t newval) { fUseUnityParticleWeights = newval; };
  void SetUseUnityEventWeights(Bool_t newval) { fUseUnityEventWeights = newval; };
  void SetSystFlag(Int_t newval) { if(!fGFWSelection) fGFWSelection = new AliGFWCuts(); fGFWSelection->SetupCuts(newval); }; //Flag for systematics
  void SetDCAxyFunctionalForm(TString newval) { fDCAxyFunctionalForm = newval; } //Call after SystFlag
  void SetConsistencyFlag(UInt_t newval) { fConsistencyFlag = newval; };
  void SetCentralityEstimator(TString newval) { if(fCentEst) delete fCentEst; fCentEst = new TString(newval); };
  void SetContSubfix(TString newval) {if(fContSubfix) delete fContSubfix; fContSubfix = new TString(newval); };
  void OverrideMCFlag(Bool_t newval) { fIsMC = newval; };
  Int_t GetNtotTracks(AliAODEvent*, const Double_t &ptmin, const Double_t &ptmax, Double_t *vtxp);
  Int_t GetNtotMCTracks(const Double_t &ptmin, const Double_t &ptmax);
  void SetUseRecoNchForMC(Bool_t newval) { fUseRecoNchForMC = newval; };
  void SetNBootstrapProfiles(Int_t newval) {if(newval<0) {printf("Number of subprofiles cannot be < 0!\n"); return; }; fNBootstrapProfiles = newval; };
  void SetPseudoEfficiency(Double_t newval) {fPseudoEfficiency = newval; };
  void SetBypassTriggerAndEventCuts(Bool_t newval) { fBypassTriggerAndEventCuts = newval; };
  void SetV0PUCut(TString newval) { if(fV0CutPU) delete fV0CutPU; fV0CutPU = new TF1("fV0CutPU", newval.Data(), 0, 100000); };
  void SetDisablePileup(bool disable) { fDisablePileup = disable; };
  void SetEnableFB768DCAxy(bool newval) { fEnableFB768dcaxy = newval;}
  void SetUseOldPileup(bool newval) { fUseOldPileup = newval; }
  void SetCentralPileup(double newval) {fCentralPU = newval;}
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskGammaSoft(const AliAnalysisTaskGammaSoft&);
  AliAnalysisTaskGammaSoft& operator=(const AliAnalysisTaskGammaSoft&);
  Int_t fSystFlag;
  Int_t fEventCutFlag; //0 for standard AliEventCuts; 1 for LHC15o pass2; 2 for LHC18qr pass3
  TString *fContSubfix;
  TString *fCentEst;
  Bool_t fIsMC;
  Bool_t fBypassTriggerAndEventCuts;
  Bool_t fDisablePileup;
  Bool_t fUseOldPileup;
  TString fDCAxyFunctionalForm;
  Bool_t fUseRecoNchForMC; //Flag to use Nch from reconstructed, when running MC closure
  TRandom *fRndm;
  Int_t fNBootstrapProfiles; //Number of profiles for bootstrapping
  TAxis *fPtAxis;
  TAxis *fEtaAxis;
  TAxis *fMultiAxis;      //Multiplicity axis (either for V0M or Nch)
  TAxis *fV0MMultiAxis;   //Defaults V0M bins
  TAxis *fDCAAxis;
  Double_t *fPtBins; //!
  Int_t fNPtBins; //!
  Double_t *fEtaBins; //!
  Int_t fNEtaBins; //!
  Double_t *fMultiBins; //!
  Int_t fNMultiBins; //!
  Double_t *fDCABins; //!
  Int_t fNDCABins; //!
  Double_t *fV0MBinsDefault; //!
  Int_t fNV0MBinsDefault; //!
  Bool_t fUseNch;
  Bool_t fUseUnityParticleWeights;
  Bool_t fUseUnityEventWeights;
  Double_t fEtaMpt;
  Double_t fEtaAcceptance;
  Double_t fEtaV2Sep;
  TList *fQAList; //
  TH1D* fEventCount; //!
  TH1D *fMultiDist;
  TH2D **fMultiVsV0MCorr; //!
  TH2D *fNchTrueVsReco; //!
  TH2D *fESDvsFB128;
  TList *fPtList;
  TList *fCovList;
  AliProfileBS *fv24deltapt2_v24mpt2; //!
  AliProfileBS *fv24deltapt2_v24mpt; //!
  AliProfileBS *fv24; //!
  AliProfileBS *fdeltapt2_mpt2; //!
  AliProfileBS *fdeltapt2_mpt; //!
  AliProfileBS *fmmpt; //!
  AliProfileBS *fv2deltapt2_v2mpt2; //!
  AliProfileBS *fv2deltapt2_v2mpt; //!
  AliProfileBS *fv2; //!
  AliProfileBS *fv2deltapt_v2mpt; //!
  UInt_t fTriggerType;
  TList *fWeightList; //!
  AliGFWWeights **fWeights;//! This should be stored in TList
  Int_t fRunNo; //!
  AliGFWCuts *fGFWSelection;
  AliGFWFlowContainer *fFC; //!
  AliGFW *fGFW; //! not stored
  vector<AliGFW::CorrConfig> corrconfigs; //! do not store
  TList *fEfficiencyList;
  vector<vector<TH2D*>> fEfficiency; //TH2Ds for efficiency calculation
  TH1D **fEfficiencies; //TH1Ds for picking up efficiencies
  Double_t fPseudoEfficiency; //Pseudo efficiency to reject tracks. Default value set to 2, only used when the value is <1
  TH3D *fDCAxyVsPt_noChi2; //!
  TH2D *fWithinDCAvsPt_withChi2; //!
  TH3D *fDCAxyVsPt_withChi2; //!
  TH2D *fWithinDCAvsPt_noChi2; //!
  TH1D *fV0MMulti; //!
  TH2I *fITSvsTPCMulti; //!
  Double_t fCorrPar[2]; //Yes need to store
  Bool_t fUseCorrCuts; //Yes need to store
  TF1 *fSPDCutPU; //Store these
  TF1 *fV0CutPU; //Store these
  TF1 *fCenCutLowPU; //Store these
  TF1 *fCenCutHighPU; //Store these
  TF1 *fMultCutPU; //Store these
  Double_t fCentralPU;
  TH3D** fPhiEtaVz; //!
  TH3D** fptDCAxyDCAz; //!
  TH1D** fChi2TPCcls; //!
  TH1D** fTPCcls; //!
  TH1D** fTPCcrsrws; //!
  TH2D* fhQAEventsfMult32vsCentr; //!
  TH2D* fhQAEventsMult128vsCentr; //!
  TH2D* fhQAEventsfMultTPCvsTOF; //!
  TH2D* fhQAEventsfMultTPCvsESD; //!
  unsigned int fEventWeight;
  vector<vector<double>>  wp;
  void DCAxyz(const AliAODTrack *track, const AliVEvent *evt, Double_t (&dcaxyz)[2]);
  Bool_t FillFCs(const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &rndmn, const Bool_t deubg=kFALSE);
  Bool_t FillCovariance(AliProfileBS* target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt, const Double_t &l_rndm);
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*, const Double_t &ptMin=0.5, const Double_t &ptMax=2, Double_t *vtxp=0);
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp, Int_t &nTot);
  void SetupAxes();
  Bool_t AcceptCustomEvent(AliAODEvent*);
  Bool_t LoadWeights(const Int_t &runno);
  UInt_t fConsistencyFlag;
  Bool_t fRequireReloadOnRunChange;
  Bool_t fEnableFB768dcaxy;
  Double_t *GetBinsFromAxis(TAxis *inax);

  ClassDef(AliAnalysisTaskGammaSoft,3);
};

#endif
