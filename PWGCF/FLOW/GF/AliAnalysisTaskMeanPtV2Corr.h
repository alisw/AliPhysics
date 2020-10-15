#ifndef MEANPTV2CORRELATIONS__H
#define MEANPTV2CORRELATIONS__H
#include "AliAnalysisTaskSE.h"
#include "TComplex.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliGFW.h"
#include "AliPID.h"

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
  AliAnalysisTaskMeanPtV2Corr(const char *name, Bool_t IsMC=kTRUE, TString StageSwitch="");
  virtual ~AliAnalysisTaskMeanPtV2Corr();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  Bool_t CheckTrigger(Double_t);
  Bool_t AcceptAOD(AliAODEvent*, Double_t lvtxXYZ[3]);
  Bool_t AcceptParticle(AliVParticle*);
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
  void FillWeights(AliAODEvent*, Double_t vz, Double_t l_Cent);
  void FillMeanPtCounter(Double_t l_pt, Double_t &l_sum, Double_t &l_count, AliGFWWeights *inWeight); //passing by ref., considering how ofter this is called
  void FillMeanPt(AliAODEvent*, Double_t vz, Double_t l_Cent);
  void FillCK(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent);
  void ProduceALICEPublished_MptProd(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent);
  void ProduceALICEPublished_CovProd(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent);
  void ProduceFBSpectra(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent);
  Int_t GetStageSwitch(TString instr);
  AliGFW::CorrConfig GetConf(TString head, TString desc, Bool_t ptdif) { return fGFW->GetCorrelatorConfig(desc,head,ptdif);};
  void CreateCorrConfigs();
  void LoadWeightAndMPT();
  void GetSingleWeightFromList(AliGFWWeights **inWeights, TString pf="");
  Bool_t WithinSigma(Double_t SigmaCut, AliAODTrack *inTrack, AliPID::EParticleType partType);
  void FillWPCounter(Double_t[5], Double_t, Double_t);
  void CalculateMptValues(Double_t[4], Double_t[5]);
  Bool_t LoadMyWeights(Int_t lRunNo = 0);
  Int_t GetBayesPIDIndex(AliAODTrack*);
  Double_t GetMyWeight(Double_t eta, Double_t phi, Int_t pidind);
  void ChangeMptSet(Bool_t newval) {fmptSet = newval; };
  void SetTrackFilterBit(Int_t newval) { fFilterBit = newval; };
  void SetDisablePID(Bool_t newval) { fDisablePID = newval; };
  void SetPtBins(Int_t nBins, Double_t *ptbins);
  void SetMultiBins(Int_t nBins, Double_t *multibins);
  void SetEta(Double_t newval) { fEta = newval; };
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskMeanPtV2Corr(const AliAnalysisTaskMeanPtV2Corr&);
  AliAnalysisTaskMeanPtV2Corr& operator=(const AliAnalysisTaskMeanPtV2Corr&);
  Int_t fStageSwitch;
  Bool_t fIsMC;
  TAxis *fPtAxis;
  TAxis *fMultiAxis;
  Double_t *fPtBins; //!
  Int_t fNPtBins; //!
  Double_t *fMultiBins; //!
  Int_t fNMultiBins; //!
  Double_t fEta;
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
  TProfile **fCovariance; //!
  Bool_t fmptSet;
  UInt_t fTriggerType; //! No need to store
  TList *fWeightList; //!
  AliGFWWeights **fWeights;//! This should be stored in TList
  TList *fNUAList; //!
  TH2D **fNUAHist; //!
  Int_t fRunNo; //!
  AliGFWCuts *fMidSelection; //!
  AliGFWCuts *fFWSelection; //!
  AliGFWFlowContainer *fFC;
  AliGFW *fGFW; //! not stored
  vector<AliGFW::CorrConfig> corrconfigs; //! do not store
  TList *fSpectraList;
  TH2D **fSpectra;
  TH1D *fV0MMulti;
  Bool_t FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndmn);
  Bool_t FillCovariance(TProfile* target, AliGFW::CorrConfig corconf, Double_t cent, Double_t d_mpt, Double_t dw_mpt);
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*,Double_t ptMin=0.5, Double_t ptMax=2, Int_t FilterBit=96);
  Int_t fFilterBit;
  Bool_t fDisablePID;
  Double_t *GetBinsFromAxis(TAxis *inax);
  ClassDef(AliAnalysisTaskMeanPtV2Corr,1);
};

#endif
