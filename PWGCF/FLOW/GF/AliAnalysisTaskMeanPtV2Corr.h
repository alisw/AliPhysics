#ifndef MEANPTV2CORRELATIONS__H
#define MEANPTV2CORRELATIONS__H
#include "AliAnalysisTaskSE.h"
#include "TComplex.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliGFW.h"
#include "AliPID.h"
#include "AliMCEvent.h"
#include "AliGFWCuts.h"
#include "TString.h"

class TList;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TProfile2D;
class TComplex;
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
  Bool_t AcceptAOD(AliAODEvent*, Double_t lvtxXYZ[3]);
  Bool_t AcceptParticle(AliVParticle*);
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
  void FillWeights(AliAODEvent*, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  void FillMeanPtCounter(Double_t l_pt, Double_t &l_sum, Double_t &l_count, AliGFWWeights *inWeight); //passing by ref., considering how ofter this is called
  void FillMeanPtCounterWW(const Double_t &l_pt, Double_t &l_sum, Double_t &l_count, const Double_t &inWeight); //passing by ref., considering how ofter this is called
  void FillMeanPt(AliAODEvent*, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  void FillMeanPtMC(AliAODEvent*, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  void FillCK(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  void ProduceALICEPublished_MptProd(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  void ProduceALICEPublished_CovProd(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  void ProduceFBSpectra(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  void ProduceEfficiencies(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  Int_t GetStageSwitch(TString instr);
  AliGFW::CorrConfig GetConf(TString head, TString desc, Bool_t ptdif) { return fGFW->GetCorrelatorConfig(desc,head,ptdif);};
  void CreateCorrConfigs();
  void LoadWeightAndMPT();
  void GetSingleWeightFromList(AliGFWWeights **inWeights, TString pf="");
  Bool_t WithinSigma(Double_t SigmaCut, AliAODTrack *inTrack, AliPID::EParticleType partType);
  void FillWPCounter(Double_t[5], Double_t, Double_t);
  void CalculateMptValues(Double_t[4], Double_t[5]);
  Bool_t LoadMyWeights(const Int_t &lRunNo = 0);
  Int_t GetBayesPIDIndex(AliAODTrack*);
  Double_t GetMyWeight(Double_t eta, Double_t phi, Int_t pidind);
  void ChangeMptSet(Bool_t newval) {fmptSet = newval; };
  Int_t GetPIDIndex(const Int_t &pdgcode);
  void SetDisablePID(Bool_t newval) { fDisablePID = newval; };
  void SetPtBins(Int_t nBins, Double_t *ptbins);
  void SetMultiBins(Int_t nBins, Double_t *multibins);
  void SetV0MBins(Int_t nBins, Double_t *multibins);
  void SetV2dPtMultiBins(Int_t nBins, Double_t *multibins);
  void SetEta(Double_t newval) { fEta = newval; };
  void SetEtaNch(Double_t newval) { fEtaNch = newval; };
  void SetEtaV2Sep(Double_t newval) { fEtaV2Sep = TMath::Abs(newval); };
  void SetUseNch(Bool_t newval) { fUseNch = newval; };
  void SetUseWeightsOne(Bool_t newval) { fUseWeightsOne = newval; };
  void ExtendV0MAcceptance(Bool_t newval) { fExtendV0MAcceptance = newval; };
  void SetSystFlag(Int_t newval) { if(!fGFWSelection) fGFWSelection = new AliGFWCuts(); fGFWSelection->SetupCuts(newval); }; //Flag for systematics
  void SetConsistencyFlag(UInt_t newval) { fConsistencyFlag = newval; };
  void SetCentralityEstimator(TString newval) { if(fCentEst) delete fCentEst; fCentEst = new TString(newval); };
  void SetContSubfix(TString newval) {if(fContSubfix) delete fContSubfix; fContSubfix = new TString(newval); if(!fContSubfix->IsNull()) fContSubfix->Prepend("_"); };
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskMeanPtV2Corr(const AliAnalysisTaskMeanPtV2Corr&);
  AliAnalysisTaskMeanPtV2Corr& operator=(const AliAnalysisTaskMeanPtV2Corr&);
  Int_t fStageSwitch;
  TString *fContSubfix;
  TString *fCentEst;
  Bool_t fExtendV0MAcceptance;
  Bool_t fIsMC;
  AliMCEvent *fMCEvent; //! MC event
  TAxis *fPtAxis;
  TAxis *fMultiAxis;
  TAxis *fV0MMultiAxis;
  Double_t *fPtBins; //!
  Int_t fNPtBins; //!
  Double_t *fMultiBins; //!
  Int_t fNMultiBins; //!
  Bool_t fUseNch;
  Bool_t fUseWeightsOne;
  Double_t fEta;
  Double_t fEtaNch;
  Double_t fEtaV2Sep; //Please don't add multiple wagons with dif. values; implement subevents in the code instead. This would save TONS of CPU time.
  AliPIDResponse *fPIDResponse; //!
  AliPIDCombined *fBayesPID; //!
  TList *fMPTList; //!
  TProfile **fmPT; //!
  TH1D *fMultiDist;
  TProfile *fNchVsMulti;
  TProfile *fNchInBins;
  TList *fptVarList;
  TProfile **fptvar; //!
  TList *fCovList;
  TList *fV2dPtList;
  TProfile **fCovariance; //!
  Bool_t fmptSet;
  UInt_t fTriggerType;
  TList *fWeightList; //!
  AliGFWWeights **fWeights;//! This should be stored in TList
  TList *fNUAList; //!
  TH2D **fNUAHist; //!
  Int_t fRunNo; //!
  AliGFWCuts *fGFWSelection;
  AliGFWFlowContainer *fFC;
  AliGFW *fGFW; //! not stored
  vector<AliGFW::CorrConfig> corrconfigs; //! do not store
  TList *fSpectraList;
  TH2D **fSpectra;
  TList *fEfficiencyList;
  TH2D **fEfficiency; //TH2Ds for efficiency calculation
  TH1D **fEfficiencies; //TH1Ds for picking up efficiencies
  TH1D *fV0MMulti;
  TH1D *fV2dPtMulti;
  Bool_t FillFCs(const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &rndmn);
  Bool_t Fillv2dPtFCs(const AliGFW::CorrConfig &corconf, const Double_t &dpt, const Double_t &rndmn, const Int_t index);
  Bool_t FillCovariance(TProfile* target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt);
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*, const Double_t &ptMin=0.5, const Double_t &ptMax=2, Double_t *vtxp=0);
  Bool_t fDisablePID;
  UInt_t fConsistencyFlag;
  Bool_t fRequireReloadOnRunChange;
  Double_t *GetBinsFromAxis(TAxis *inax);
  ClassDef(AliAnalysisTaskMeanPtV2Corr,1);
};

#endif
