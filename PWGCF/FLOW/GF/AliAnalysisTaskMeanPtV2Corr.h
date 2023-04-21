//Author: Vytautas Vislavicius
#ifndef MEANPTV2CORRELATIONS__H
#define MEANPTV2CORRELATIONS__H
#include "AliAnalysisTaskSE.h"
// #include "TComplex.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliGFW.h"
#include "AliPID.h"
#include "AliMCEvent.h"
#include "AliGFWCuts.h"
#include "TString.h"
#include "AliProfileBS.h"
#include "AliCkContainer.h"
#include "TRandom.h"
#include "AliMCSpectraWeights.h"
#include "AliAODTracklets.h"
#include "AliAODVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliVMultiplicity.h"
#include "AliGFWFilter.h"
#include "TRegexp.h"
#include "TFile.h"


class TList;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TProfile2D;
// class TComplex;
class AliAODEvent;
class AliVTrack;
class AliVVertex;
class AliInputEventHandler;
class AliAODTrack;
class TClonesArray;
class AliAODVertex;
class AliAnalysisUtils;
class TProfile;
class AliGFWWeights;
class AliVParticle;
class AliGFWCuts;
class AliGFWFlowContainer;
class AliPIDResponse;
class AliPIDCombined;
class AliEffFDContainer;

class AliAnalysisTaskMeanPtV2Corr : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMeanPtV2Corr();
  AliAnalysisTaskMeanPtV2Corr(const char *name, Bool_t IsMC=kTRUE, TString StageSwitch="", TString ContainerSubfix="");
  virtual ~AliAnalysisTaskMeanPtV2Corr();
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  Bool_t CheckTrigger(Double_t);
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
  void SetEventCutFlag(Int_t newval) { fEventCutFlag = newval; };
  void SetupFlagsByIndex(Int_t); //Setting up event and track flags
  void CovSkipMpt(AliGFWFlags *lFlags, AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  Int_t GetStageSwitch(TString instr);
  AliGFW::CorrConfig GetConf(TString head, TString desc, Bool_t ptdif) { return fGFW->GetCorrelatorConfig(desc.Data(),head.Data(),ptdif);};
  void CreateCorrConfigs();
  void FillWPCounter(Double_t[5], Double_t, Double_t);
  Bool_t LoadMyWeights(const Int_t &lRunNo = 0);
  Int_t GetBayesPIDIndex(AliVTrack*);
  Int_t GetPIDIndex(const Int_t &pdgcode);
  void SetDisablePID(Bool_t newval) { fDisablePID = newval; };
  void SetPtBins(Int_t nBins, Double_t *ptbins);
  void SetMultiBins(Int_t nBins, Double_t *multibins);
  void SetV0MBins(Int_t nBins, Double_t *multibins);
  void SetEtaEffBins(Int_t nBins, Double_t *etabins);
  void SetV2dPtMultiBins(Int_t nBins, Double_t *multibins);
  void SetEta(Double_t newval) { fEta = newval; fEtaLow=-9999; };
  void SetEta(Double_t etaLow, Double_t etaHigh) { fEtaLow = etaLow; fEta = etaHigh; };
  void SetEtaNch(Double_t newval) { fEtaNch = newval; };
  void SetEtaV2Sep(Double_t newval) { fEtaV2Sep = newval; };
  void SetUseNch(Bool_t newval) { fUseNch = newval; };
  void SetUseWeightsOne(Bool_t newval) { fUseWeightsOne = newval; };
  void ExtendV0MAcceptance(Bool_t newval) { fExtendV0MAcceptance = newval; };
  void SetSystFlag(Int_t newval) { if(!fGFWSelection) fGFWSelection = new AliGFWCuts(); fGFWSelection->SetupCuts(newval); }; //Flag for systematics
  void SetConsistencyFlag(UInt_t newval) { fConsistencyFlag = newval; };
  void SetCentralityEstimator(TString newval) { if(fCentEst) delete fCentEst; fCentEst = new TString(newval); };
  void SetContSubfix(TString newval) {if(fContSubfix) delete fContSubfix; fContSubfix = new TString(newval); };
  void OverrideMCFlag(Bool_t newval) { fIsMC = newval; };
  Int_t GetNtotTracks(AliAODEvent*, const Double_t &ptmin, const Double_t &ptmax, Double_t *vtxp);
  void SetUseRecoNchForMC(Bool_t newval) { fUseRecoNchForMC = newval; };
  void SetNBootstrapProfiles(Int_t newval) {if(newval<0) {printf("Number of subprofiles cannot be < 0!\n"); return; }; fNBootstrapProfiles = newval; };
  void SetWeightSubfix(TString newval) { fWeightSubfix=newval; }; //base (runno) + subfix (systflag), delimited by ;. First argument always base, unless is blank. In that case, w{RunNo} is used for base.
  void SetPseudoEfficiency(Double_t newval) {fPseudoEfficiency = newval; };
  void SetBypassTriggerAndEventCuts(Bool_t newval) { fBypassTriggerAndEvetCuts = newval; };
  void SetV0PUCut(TString newval) { if(fV0CutPU) delete fV0CutPU; fV0CutPU = new TF1("fV0CutPU", newval.Data(), 0, 100000);
}
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskMeanPtV2Corr(const AliAnalysisTaskMeanPtV2Corr&);
  AliAnalysisTaskMeanPtV2Corr& operator=(const AliAnalysisTaskMeanPtV2Corr&);
  Int_t fStageSwitch;
  Int_t fSystFlag;
  Int_t fEventCutFlag; //0 for standard AliEventCuts; 1 for LHC15o pass2; 2 for LHC18qr pass3
  Int_t fEvNomFlag; //GFWFilter implementation
  Int_t fTrNomFlag; //GFWFilter implementation
  TString *fContSubfix;
  TString *fCentEst;
  Bool_t fExtendV0MAcceptance;
  Bool_t fIsMC;
  Bool_t fBypassTriggerAndEvetCuts;
  AliMCEvent *fMCEvent; //! MC event
  Bool_t fUseRecoNchForMC; //Flag to use Nch from reconstructed, when running MC closure
  TRandom *fRndm; //For random number generation
  Int_t fNBootstrapProfiles; //Number of profiles for bootstrapping
  TAxis *fPtAxis;
  TAxis *fMultiAxis;
  TAxis *fV0MMultiAxis;
  TAxis *fEtaAxis; //Only relevant for efficiencies when running p-Pb
  Double_t *fPtBins; //!
  Int_t fNPtBins; //!
  Double_t *fMultiBins; //!
  Int_t fNMultiBins; //!
  Double_t *fEtaBins; //!
  Int_t fNEtaBins; //!
  Bool_t fUseNch;
  Bool_t fUseWeightsOne;
  Double_t fEta;
  Double_t fEtaLow;
  Double_t fEtaNch;
  Double_t fEtaV2Sep; //Please don't add multiple wagons with dif. values; implement subevents in the code instead. This would save TONS of CPU time.
  AliPIDResponse *fPIDResponse; //!
  AliPIDCombined *fBayesPID; //!
  TList *fQAList; //
  TH1D *fMultiDist;
  TH2D **fMultiVsV0MCorr; //!
  TH2D *fNchTrueVsReco; //!
  TH2D *fESDvsFB128;
  TProfile *fNchVsMulti;
  TProfile *fNchInBins;
  TList *fptVarList;
  AliCkContainer *fCkCont;
  TList *fCovList;
  TList *fV2dPtList;
  AliProfileBS **fCovariance; //!
  UInt_t fTriggerType;
  TList *fWeightList; //!
  AliGFWWeights **fWeights;//! This should be stored in TList
  TString fWeightSubfix;
  Int_t fRunNo; //!
  AliGFWCuts *fGFWSelection;
  AliGFWCuts *fGFWNtotSelection;
  AliGFWFlowContainer *fFC;
  AliGFW *fGFW; //! not stored
  vector<AliGFW::CorrConfig> corrconfigs; //! do not store
  TList *fEfficiencyList;
  TH2D **fEfficiency; //TH2Ds for efficiency calculation
  TH1D **fEfficiencies; //TH1Ds for picking up efficiencies
  Double_t fPseudoEfficiency; //Pseudo efficiency to reject tracks. Default value set to 2, only used when the value is <1
  TH3D *fDCAxyVsPt_noChi2;
  TH2D *fWithinDCAvsPt_withChi2;
  TH3D *fDCAxyVsPt_withChi2;
  TH2D *fWithinDCAvsPt_noChi2;
  TH1D *fV0MMulti;
  TH2D *fITSvsTPCMulti;
  TH1D *fV2dPtMulti;
  Double_t fCorrPar[2]; //Yes need to store
  Bool_t fUseCorrCuts; //Yes need to store
  TF1 *fSPDCutPU; //Store these
  TF1 *fV0CutPU; //Store these
  TF1 *fCenCutLowPU; //Store these
  TF1 *fCenCutHighPU; //Store these
  TF1 *fMultCutPU; //Store these
  AliESDtrackCuts *fStdTPCITS2011; //Needed for counting tracks for custom event cuts
  Bool_t LoadWeights(const Int_t &runno);
  Bool_t FillFCs(const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &rndmn, const Bool_t deubg=kFALSE);
  Bool_t Fillv2dPtFCs(const AliGFW::CorrConfig &corconf, const Double_t &dpt, const Double_t &rndmn, const Int_t index);
  Bool_t FillCovariance(AliProfileBS* target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt, const Double_t &l_rndm);
  Bool_t fDisablePID;
  UInt_t fConsistencyFlag;
  Bool_t fRequireReloadOnRunChange;
  Double_t *GetBinsFromAxis(TAxis *inax);
  ClassDef(AliAnalysisTaskMeanPtV2Corr,1);
};

#endif
