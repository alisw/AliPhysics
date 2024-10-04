#ifndef PTCORR__H
#define PTCORR__H

#include "AliAnalysisTaskSE.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
#include "AliEventCuts.h"
#include "TString.h"
#include "TRandom.h"
#include "TAxis.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "AliPtPtContainer.h"
#include "AliProfileBS.h"
//#include "AliProfile2DBS.h"
#include "TProfile2D.h"
#include "TH3D.h"
#include "AliGFWCuts.h"
#include "TF1.h"
#include "AliESDtrackCuts.h"

class AliAODEvent;
class AliAODTrack;
class AliESDEvent;


class AliAnalysisTaskPtCorr : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPtCorr();
  AliAnalysisTaskPtCorr(const char *name, Bool_t IsMC=kTRUE, TString ContainerSubfix="");
  virtual ~AliAnalysisTaskPtCorr();
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  AliMCEvent *getMCEvent();
  Bool_t CheckTrigger(Double_t);
  Bool_t AcceptAOD(AliAODEvent*, Double_t lvtxXYZ[3]);
  Bool_t IsPileupEvent(AliAODEvent* ev, double centrality);
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
  void LoadCorrectionsFromLists();
  void ProcessOnTheFly();
  void ProcessGen();
  Bool_t IsMCEventSelected(TObject *obj);
  void FillWPCounter(vector<vector<double>> &inarr, double w, double p);
  void SetFillCM(unsigned int flag) { fCMflag = flag; }
  void SetPtBins(Int_t nBins, Double_t *ptbins);
  void SetEtaBins(Int_t nBins, Double_t *etabins);
  void SetMultiBins(Int_t nBins, Double_t *multibins);
  void SetMultiBins(Int_t nBins, Double_t low, Double_t high);
  void SetV0MBins(Int_t nBins, Double_t *multibins);
  void SetMptBins(Int_t nMptBins, Double_t *mptbins);
  void SetMptBins(Int_t nMptBins, Double_t mptlow, Double_t mpthigh);
  void SetDCABins(Int_t nBins, Double_t *dcabins);
  void SetDCABins(Int_t nBins, Double_t low, Double_t high);
  void SetEtaAcceptance(Double_t eta) { fEtaMptAcceptance[0] = -TMath::Abs(eta); fEtaMptAcceptance[1] = TMath::Abs(eta); fEtaMultAcceptance[0] = -TMath::Abs(eta); fEtaMultAcceptance[1] = TMath::Abs(eta);}
  void SetEtaMptAcceptance(Double_t low, Double_t up) { fEtaMptAcceptance[0] = low; fEtaMptAcceptance[1] = up; };
  void SetEtaMultAcceptance(Double_t low, Double_t up) { fEtaMultAcceptance[0] = low; fEtaMultAcceptance[1] = up; };
  void SetEtaAbsolute(Bool_t newval) { fEtaAbsolute = newval; }
  void SetUseNch(Bool_t newval) { fUseNch = newval; };
  void SetUseWeightsOne(Bool_t newvalNUE) { fUseNUEOne = newvalNUE; };
  void SetSystFlag(Int_t newval) { if(!fGFWSelection) fGFWSelection = new AliGFWCuts(); fGFWSelection->SetupCuts(newval); }; //Flag for systematics
  void SetDCAxyFunctionalForm(TString newval) { fDCAxyFunctionalForm = newval; } //Call after SystFlag
  void SetCentralityEstimator(TString newval) { if(fCentEst) delete fCentEst; fCentEst = new TString(newval); };
  void SetContSubfix(TString newval) {if(fContSubfix) delete fContSubfix; fContSubfix = new TString(newval); };
  void OverrideMCFlag(Bool_t newval) { fIsMC = newval; };
  int GetNtotTracks(AliAODEvent*, const Double_t &ptmin, const Double_t &ptmax, Double_t *vtxp, Int_t iCent);
  int GetNtotMCTracks(const Double_t &ptmin, const Double_t &ptmax);
  void SetUseRecoNchForMC(Bool_t newval) { fUseRecoNchForMC = newval; };
  void SetNBootstrapProfiles(Int_t newval) {if(newval<0) {printf("Number of subprofiles cannot be < 0!\n"); return; }; fNBootstrapProfiles = newval; };
  void SetWeightSubfix(TString newval) { fWeightSubfix=newval; }; //base (runno) + subfix (systflag), delimited by ;. First argument always base, unless is blank. In that case, w{RunNo} is used for base.
  void SetPseudoEfficiency(Double_t newval) {fPseudoEfficiency = newval; };
  void SetBypassTriggerAndEventCuts(Bool_t newval) { fBypassTriggerAndEventCuts = newval; };
  void SetV0PUCut(TString newval) { if(fV0CutPU) delete fV0CutPU; fV0CutPU = new TF1("fV0CutPU", newval.Data(), 0, 100000); };
  void SetEventWeight(unsigned int weight) { fEventWeight = weight; };
  void SetDisablePileup(bool disable) { fDisablePileup = disable; };
  void SetFillQA(Bool_t newval) { fFillQA = newval; }
  void SetPtMPar(int newval) { fPtMpar = newval; }
  void SetEnableFB768DCAxy(bool newval) { fEnableFB768dcaxy = newval;}
  void SetUseOldPileup(bool newval) { fUseOldPileup = newval; }
  void SetCentralPileup(double newval) {fCentralPU = newval;}
  void SetOnTheFly(bool newval) {fOnTheFly = newval;}
  void SetOnTheFlyGen(bool newval) {fOnTheFlyGen = newval;}
  double getGeneratorCentrality();
  void SetOTFGenerator(TString gen) { fGenerator = gen; }
  void SetUseIP(bool newval) { fUseIP = newval;}
  void SetUseCentCalibration(bool newval) { fUseCentCalibration = newval; }
  void SetRejectMBTriggeredEventsMarkedSpecialTrigger(bool newval) { fRejectMBtriggerEventsMarkedSpecial = newval; }
  void SwitchV0MToETmult(bool newval) { fUseEToverV0M = newval; }
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskPtCorr(const AliAnalysisTaskPtCorr&);
  AliAnalysisTaskPtCorr& operator=(const AliAnalysisTaskPtCorr&);
  Int_t fSystFlag;
  Int_t fEventCutFlag; //0 for standard AliEventCuts; 1 for LHC15o pass2; 2 for LHC18qr pass3
  TString *fContSubfix;
  TString *fCentEst;
  Bool_t fIsMC;
  Bool_t fBypassTriggerAndEventCuts;
  Bool_t fDisablePileup;
  Bool_t fUseOldPileup;
  Bool_t fCentSelectForMptNch;
  Bool_t fUseIP;
  Bool_t fUseCentCalibration;
  Bool_t fRejectMBtriggerEventsMarkedSpecial;
  Bool_t fUseEToverV0M;
  unsigned int fCMflag;
  TString fDCAxyFunctionalForm;
  Bool_t fOnTheFly;
  Bool_t fOnTheFlyGen;
  TString fGenerator;
  AliMCEvent *fMCEvent; //!<! MC event
  AliInputEventHandler* fMCHandler;  //!<! MCEventHandler
  AliStack* fStack;
  Bool_t fUseRecoNchForMC; //Flag to use Nch from reconstructed, when running MC closure
  TRandom *fRndm;
  Int_t fNBootstrapProfiles; //Number of profiles for bootstrapping
  Bool_t fFillQA;
  TAxis *fPtAxis;
  TAxis *fDCAAxis;
  TAxis *fEtaAxis;
  TAxis *fMultiAxis;      //Multiplicity axis (either for V0M or Nch)
  TAxis *fCentPtAxis;      //Centrality axis for mpt fluctuations
  TAxis *fMptAxis;      //Mpt axis
  TAxis *fV0MMultiAxis;   //Defaults V0M bins
  Double_t *fPtBins; //!
  Int_t fNPtBins; //!
  Double_t *fDCABins; //!
  Int_t fNDCABins; //!
  Double_t *fEtaBins; //!
  Int_t fNEtaBins; //!
  Double_t *fMultiBins; //!
  Int_t fNMultiBins; //!
  Double_t *fCentPtBins; //!
  Int_t fNCentPtBins; //!
  Double_t *fMptBins; //!
  Int_t fNMptBins; //!
  Double_t *fV0MBinsDefault; //!
  Int_t fNV0MBinsDefault; //!
  Bool_t fUseNch;
  Bool_t fUseNUEOne;
  Int_t fPtMpar;
  Double_t fEtaLow;
  vector<Double_t> fEtaMptAcceptance;
  vector<Double_t> fEtaMultAcceptance;
  Bool_t fEtaAbsolute;
  Double_t fImpactParameterMC;
  TList *fQAList; //
  TH1D* fEventCount; //!
  TH1D* fV0MMulti; //!
  TH2D *fNchTrueVsReco; //!
  TH2D *fESDvsFB128;
  TList *fptList;
  AliPtPtContainer  *fPtCont;
  AliPtPtContainer  *fPtContV0Mmult;
  AliPtPtContainer  *fPtContET;
  AliPtPtContainer  *fPtContCent;
  UInt_t fTriggerType;
  TList *fSpectraList; //!
  TH3D **fSpectraGen; //!
  TH3D **fSpectraRec; //!
  TH2D **fDetectorResponse; //!
  TString fWeightSubfix;
  Int_t fRunNo; //!
  AliGFWCuts *fGFWSelection;
  AliGFWCuts *fGFWNtotSelection;
  TList *fEfficiencyList;
  vector<vector<TH2D*>> fEfficiency; //TH2Ds for efficiency calculation
  TH1D **fEfficiencies; //TH1Ds for picking up efficiencies
  TH1* fCentcal; //TH1 for OTF centrality calibration
  Double_t fPseudoEfficiency; //Pseudo efficiency to reject tracks. Default value set to 2, only used when the value is <1
  TH1D* hMulti; //!
  TH2D* hCentVsNch; //!
  TH2D* hETvsNch; //!
  TH2D* hV0MvsNch; //!
  TProfile* pMeanNchVsNch; //!
  TProfile* pMeanNchVsCent; //!
  TProfile* pMeanNchvsET; //!
  TProfile* pMeanNchvsV0M; //!
  TH3D *fptDCAxyDCAz; //!
  TH3D *fPhiEtaVtxZ; //!
  TH1D* fIP; //!
  Double_t fCorrPar[2]; //Yes need to store
  Bool_t fUseCorrCuts; //Yes need to store
  TF1 *fSPDCutPU; //Store these
  TF1 *fV0CutPU; //Store these
  TF1 *fCenCutLowPU; //Store these
  TF1 *fCenCutHighPU; //Store these
  TF1 *fMultCutPU; //Store these
  Double_t fCentralPU;
  TH2D* fhQAEventsfMult32vsCentr; //!
  TH2D* fhQAEventsMult128vsCentr; //!
  TH2D* fhQAEventsfMultTPCvsTOF; //!
  TH2D* fhQAEventsfMultTPCvsESD; //!
  TH1D* fEtaMpt; //!
  TH1D* fEtaMult; //!
  int EventNo;
  unsigned int fEventWeight;
  vector<vector<double>>  wp;
  vector<vector<double>>  wpcm;
  std::map<double,double> centralitymap;
  AliESDtrackCuts *fStdTPCITS2011; //Needed for counting tracks for custom event cuts
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*, const Double_t &ptMin=0.5, const Double_t &ptMax=2, Double_t *vtxp=0);
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp, Int_t iCent, Double_t &nTot);
  void DCAxyz(const AliAODTrack *track, const AliVEvent *evt,Double_t (&dcaxyz)[2]);
  void FillAdditionalTrackQAPlots(AliAODTrack &track, const Double_t &cent, Double_t weff, Double_t wacc, const Double_t &vz, Double_t* vtxp, Bool_t beforeCuts);
  void SetupAxes();
  Bool_t AcceptCustomEvent(AliAODEvent*);
  Bool_t AcceptCustomEvent(AliESDEvent*);
  Bool_t fRequireReloadOnRunChange;
  Bool_t fEnableFB768dcaxy;
  Double_t *GetBinsFromAxis(TAxis *inax);

  ClassDef(AliAnalysisTaskPtCorr,4);
};

#endif
