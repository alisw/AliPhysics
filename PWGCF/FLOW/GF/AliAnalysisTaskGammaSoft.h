#ifndef PWGCF_FLOW_GF_GAMMASOFT_H_
#define PWGCF_FLOW_GF_GAMMASOFT_H_
#include "AliAnalysisTaskSE.h"
#include "TComplex.h"
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
#include "AliAODTracklets.h"
#include "AliAODVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliVMultiplicity.h"
#include "AliPtPtContainer.h"


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

class AliAnalysisTaskGammaSoft : public AliAnalysisTaskSE {
 public:
  enum GENERATOR {kAMPT,kHIJING};
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
  //void ProduceEfficiencies(AliESDEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  void ProcessTracks(AliAODEvent *fAOD, const Double_t &vz, const Double_t &l_Cent, Double_t *vtxp);
  AliGFW::CorrConfig GetConf(TString head, TString desc, Bool_t ptdif) { return fGFW->GetCorrelatorConfig(desc.Data(),head.Data(),ptdif);};
  void CreateCorrConfigs();
  void LoadCorrectionsFromLists();
  void FillWPCounter(vector<vector<double>> &inarr, double w, double p);
  void SetPtBins(Int_t nBins, Double_t *ptbins);
  void SetEtaBins(Int_t nBins, Double_t *etabins);
  void SetMultiBins(Int_t nBins, Double_t *multibins);
  void SetV0MBins(Int_t nBins, Double_t *multibins);
  void SetEtaMpt(Double_t newval) { fEtaMpt = newval; };
  void SetEtaAcceptance(Double_t newval) { fEtaAcceptance = newval; };
  void SetEtaV2Sep(Double_t newval) { fEtaV2Sep = newval; };
  void SetUseNch(Bool_t newval) { fUseNch = newval; };
  void SetUseWeightsOne(Bool_t newvalNUA, Bool_t newvalNUE) { fUseNUAOne = newvalNUA; fUseNUEOne = newvalNUE; };
  void SetUseEventWeightsOne(Bool_t newval) { fUseEventWeightOne = newval; };
  void SetUseStandardMethod(bool newval) { fFillStdMethod = newval; }
  void ExtendV0MAcceptance(Bool_t newval) { fExtendV0MAcceptance = newval; };
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
  void SetWeightSubfix(TString newval) { fWeightSubfix=newval; }; //base (runno) + subfix (systflag), delimited by ;. First argument always base, unless is blank. In that case, w{RunNo} is used for base.
  void SetPseudoEfficiency(Double_t newval) {fPseudoEfficiency = newval; };
  Bool_t CheckNchCorrelation(const Int_t &lNchGen, const Int_t &lNchRec) { return (fCorrPar[0]*lNchGen + fCorrPar[1] < lNchRec); };
  void SetBypassTriggerAndEventCuts(Bool_t newval) { fBypassTriggerAndEventCuts = newval; };
  void SetV0PUCut(TString newval) { if(fV0CutPU) delete fV0CutPU; fV0CutPU = new TF1("fV0CutPU", newval.Data(), 0, 100000); };
  void SetEventWeight(unsigned int weight) { fEventWeight = weight; };
  void SetDisablePileup(bool disable) { fDisablePileup = disable; };
  void SetOnTheFly(bool newval) {fOnTheFly = newval;}
  void SetOTFGenerator(TString gen) { fGenerator = gen; }
  void SetIPBins(Int_t nBins, Double_t *multibins);
  void SetAMPTCentralityMap(vector<double> b, vector<double> cent) { for(size_t i(0); i<b.size(); ++i) centralitymap[b[i]]=cent[i]; }
  void SetFillAdditionalTrackQAPlots(Bool_t newval) { fFillAdditionalQA = newval; }
  void SetPtMPar(int newval) { fPtMpar = newval; }
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
  Bool_t fExtendV0MAcceptance;
  Bool_t fIsMC;
  Bool_t fBypassTriggerAndEventCuts;
  Bool_t fDisablePileup;
  Bool_t fUseOldPileup;
  Bool_t fFillStdMethod;
  TString fDCAxyFunctionalForm;
  Bool_t fOnTheFly;
  TString fGenerator;
  AliMCEvent *fMCEvent; //! MC event
  Bool_t fUseRecoNchForMC; //Flag to use Nch from reconstructed, when running MC closure
  TRandom *fRndm;
  Int_t fNBootstrapProfiles; //Number of profiles for bootstrapping
  Bool_t fFillAdditionalQA;
  TAxis *fPtAxis;
  TAxis *fEtaAxis;
  TAxis *fMultiAxis;      //Multiplicity axis (either for V0M or Nch)
  TAxis *fV0MMultiAxis;   //Defaults V0M bins
  Double_t *fPtBins; //!
  Int_t fNPtBins; //!
  Double_t *fEtaBins; //!
  Int_t fNEtaBins; //!
  Double_t *fMultiBins; //!
  Int_t fNMultiBins; //!
  Double_t *fV0MBinsDefault; //!
  Int_t fNV0MBinsDefault; //!
  Bool_t fUseNch;
  Bool_t fUseNUAOne;
  Bool_t fUseNUEOne;
  Bool_t fUseEventWeightOne;
  Int_t fPtMpar;
  Double_t fEtaMpt;
  Double_t fEtaLow;
  Double_t fEtaAcceptance;
  Double_t fEtaV2Sep;
  TList *fQAList; //
  TH1D* fEventCount; //!
  TH1D *fMultiDist;
  TH2D *fNchTrueVsReco; //!
  AliPtPtContainer  *fPtCont;
  TList *fCovList;
  TList *fptList;
  AliProfileBS **fCovariance; //!
  UInt_t fTriggerType;
  TList *fWeightList; //!
  AliGFWWeights **fWeights;//! This should be stored in TList
  TString fWeightSubfix;
  AliGFWCuts *fGFWSelection;
  AliGFWCuts *fGFWNtotSelection;
  AliGFWFlowContainer *fFC;
  AliGFW *fGFW; //! not stored
  vector<AliGFW::CorrConfig> corrconfigs; //! do not store
  TList *fEfficiencyList;
  vector<vector<TH2D*>> fEfficiency; //TH2Ds for efficiency calculation
  TH1D **fEfficiencies; //TH1Ds for picking up efficiencies
  Double_t fPseudoEfficiency; //Pseudo efficiency to reject tracks. Default value set to 2, only used when the value is <1
  TH3D *fDCAxyVsPt_noChi2;
  TH2D *fWithinDCAvsPt_withChi2;
  TH3D *fDCAxyVsPt_withChi2;
  TH2D *fWithinDCAvsPt_noChi2;
  TH1D *fV0MMulti;
  TH1D* fIP;
  Double_t fCorrPar[2]; //Yes need to store
  Bool_t fUseCorrCuts; //Yes need to store
  TF1 *fSPDCutPU; //Store these
  TF1 *fV0CutPU; //Store these
  TF1 *fCenCutLowPU; //Store these
  TF1 *fCenCutHighPU; //Store these
  TF1 *fMultCutPU; //Store these
  Double_t fCentralPU;
  TH3D** fPhiEtaVz; //!
  TH2D** fPt; //!
  TH2D** fDCAxy; //!
  TH2D** fDCAz; //!
  TH1D** fChi2TPCcls; //!
  TH1D** fTPCcls; //!
  TH1D* fEtaMptAcceptance; //!
  TH1D* fPtMptAcceptance; //!
  TH1D* fAcceptedNch; //!
  TH2D* fhQAEventsfMult32vsCentr; //!
  TH2D* fhQAEventsMult128vsCentr; //!
  TH2D* fhQAEventsfMultTPCvsTOF; //!
  TH2D* fhQAEventsfMultTPCvsESD; //!
  Double_t fImpactParameterMC;
  unsigned int fEventWeight;
  vector<vector<double>>  wp;
  vector<vector<vector<vector<std::complex<double>>>>> abcd;
  vector<vector<vector<vector<std::complex<double>>>>> wabcd;
  std::map<double,double> centralitymap;
  AliESDtrackCuts *fStdTPCITS2011; //Needed for counting tracks for custom event cuts
  Bool_t FillFCs(const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &rndmn, const Bool_t deubg=kFALSE);
  Bool_t FillCovariance(AliProfileBS* target, const AliGFW::CorrConfig &corconf, const Double_t &cent, const Double_t &d_mpt, const Double_t &dw_mpt, const Double_t &l_rndm);
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*, const Double_t &ptMin=0.5, const Double_t &ptMax=2, Double_t *vtxp=0);
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp, Int_t &nTot);
  void FillAdditionalTrackQAPlots(AliAODTrack &track, const Double_t &cent, Double_t weff, Double_t wacc, const Double_t &vz, Double_t* vtxp, Bool_t beforeCuts);
  template<typename T>
  void FillABCDCounter(vector<vector<vector<vector<T>>>> &inarr, T a, T b, T c, T d);
  double getStdAABBCC(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  double getStdAABBCD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  double getStdAABBDD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  double getStdAABBC(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  double getStdABCC(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  double getStdABCD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  double getStdABDD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  double getStdABC(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  double getStdABD(vector<vector<vector<vector<std::complex<double>>>>> &abcdvec);
  template<typename... args>
  std::complex<double> Q(double w, double nphi, args... wnphi);
  std::complex<double> Q(double w, double nphi);
  template <typename... args>
  double P(double w, double pt, args... wpt);
  double P(double w, double pt);
  void SetupAxes();
  void CreateVnMptOutputObjects();
  Bool_t AcceptCustomEvent(AliAODEvent*);
  Bool_t LoadWeights(const Int_t &runno);
  AliMCEvent *getMCEvent();
  double getGeneratorCentrality();
  void ProcessOnTheFly();
  UInt_t fConsistencyFlag;
  Bool_t fRequireReloadOnRunChange;
  Bool_t fEnableFB768dcaxy;
  Double_t *GetBinsFromAxis(TAxis *inax);

  ClassDef(AliAnalysisTaskGammaSoft,1);
};

#endif
