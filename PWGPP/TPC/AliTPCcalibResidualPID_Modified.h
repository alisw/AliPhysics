#ifndef AliTPCcalibResidualPID_Modified_H
#define AliTPCcalibResidualPID_Modified_H

#include "AliAnalysisTaskSE.h"
#include "TTreeStream.h"
#include "TString.h"
#include "AliInputEventHandler.h"

class TArrayF;
template <class X>
class THnSparseT;
typedef class THnSparseT<TArrayF> THnSparseF;
class TFile;
class TGraphErrors;
class AliESDEvent;
class AliMCEvent;
class AliESDtrackCuts;
class AliPIDResponse;
class AliESD;
class AliAnalysisTask;
class AliESDInputHandler;
class AliESDv0KineCuts;
class AliAnalysisManager;
class AliCentrality;
class AliAnalysisUtils;
class TTree;
class TSystem;
class TStyle;
class TROOT;
class Riostream;
class TChain;
class TH2;
class TF1;
class TH1;
class TObjArray;
class TCanvas;
class TGraph;


class AliTPCcalibResidualPID_Modified : public AliAnalysisTaskSE {
 public:
  enum FitType { kAleph = 0, kLund = 1, kSaturatedLund = 2, kAlephWithAdditionalParam = 3, kAlephExternal=4 };
  enum kParticle { kElectron = 0, kPion, kKaon, kProton };
  AliTPCcalibResidualPID_Modified();
  AliTPCcalibResidualPID_Modified(const char *name);
  virtual ~AliTPCcalibResidualPID_Modified();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0);
  virtual void   Terminate(const Option_t*);
  Int_t          CompareFloat(Float_t f1=1, Float_t f2=0) const;
  //setter
  virtual void   SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  virtual void   SetESDtrackCutsV0(AliESDtrackCuts * trackCutsV0){fESDtrackCutsV0 = trackCutsV0;};
  virtual void   SetProduceTPCsignalTHnSparse(Int_t producetpcsignal){fProduceTPCSignalSparse = producetpcsignal;};
  virtual void   SetProducePIDqa(Int_t produceGlobal){fProduceGlobal = produceGlobal;};
  virtual void   SetProduceAllPadsPID(Int_t produceAllpadTypes){fProduceAllPadTypes = produceAllpadTypes;};
  virtual void   SetProduceShortPadsPID(Int_t produceShortpads){fProduceShortPads = produceShortpads;};
  virtual void   SetProduceMediumPadsPID(Int_t produceMediumpads){fProduceMediumPads = produceMediumpads;};
  virtual void   SetProduceLongPadsPID(Int_t produceLongpads){fProduceLongPads = produceLongpads;};
  virtual void   SetProduceOrocPID(Int_t produceOroc){fProduceOroc = produceOroc;};
  
  virtual Bool_t GetVertexIsOk(AliVEvent* event) const;
  
  virtual Bool_t GetUseTPCCutMIGeo() const { return fUseTPCCutMIGeo; };
  virtual void SetUseTPCCutMIGeo(Bool_t newValue) { fUseTPCCutMIGeo = newValue; };
  
  virtual Bool_t GetIsPbpOrpPb() const { return fIsPbpOrpPb; };
  virtual void SetIsPbpOrpPb(Bool_t newValue) { fIsPbpOrpPb = newValue; };
  
  virtual Bool_t GetIsPbPb() const { return fIsPbPb; };
  virtual void SetIsPbPb(Bool_t newValue) { fIsPbPb = newValue; };
  
  Double_t GetZvtxCutEvent() const { return fZvtxCutEvent; };
  virtual void SetZvtxCutEvent(Double_t newValue) { fZvtxCutEvent = newValue; };
  
  Bool_t GetCorrectdEdxEtaDependence() const { return fCorrectdEdxEtaDependence; };
  virtual void SetCorrectdEdxEtaDependence(Bool_t flag) { fCorrectdEdxEtaDependence = flag; };

  Bool_t GetCorrectdEdxMultiplicityDependence() const { return fCorrectdEdxMultiplicityDependence; };
  virtual void SetCorrectdEdxMultiplicityDependence(Bool_t flag) { fCorrectdEdxMultiplicityDependence = flag; };

  Bool_t GetCorrectdEdxPileupDependence() const { return fCorrectdEdxPileupDependence; };
  virtual void SetCorrectdEdxPileupDependence(Bool_t flag) { fCorrectdEdxPileupDependence = flag; };

  Bool_t GetCutOnProdRadiusForV0el() const { return fCutOnProdRadiusForV0el; };
  virtual void SetCutOnProdRadiusForV0el(Bool_t flag) { fCutOnProdRadiusForV0el = flag; };
  
  virtual Char_t GetV0tag(Int_t trackIndex) const;

  Bool_t GetUseMCinfo() const { return fUseMCinfo; };
  virtual void SetUseMCinfo(Bool_t flag) { fUseMCinfo = flag; };
  
  Bool_t GetWriteAdditionalOutput() const { return fWriteAdditionalOutput; };
  virtual void   SetWriteAdditionalOutput(Bool_t flag = kTRUE) { fWriteAdditionalOutput = flag; };

  virtual Int_t GetV0motherIndex(Int_t trackIndex) const;
  virtual Int_t GetV0motherPDG(Int_t trackIndex) const;
  //
  // static functions for postprocessing
  //
  Bool_t ProcessV0Tree(TTree* tree, THnSparseF* h, const Int_t recoPass=4, const TString runList="", const Bool_t excludeRuns=kFALSE);
  THnSparseF* ProcessV0TreeFile(TString filePathName, const Int_t recoPass=4, const TString runList="", const Bool_t excludeRuns=kFALSE);

  static void CreatePlotWithOwnParameters(THnSparseF * histPidQA, const Bool_t useV0s, const Char_t * type, const Char_t * system, const Double_t* parameters, AliTPCcalibResidualPID_Modified::FitType fitType, Float_t from = 0.9, Float_t to = 10e4);
  static Double_t* ExtractResidualPID(THnSparseF * histPidQA,
                                      const Bool_t useV0s = kTRUE,
                                      const Char_t * outFile = "out.root",
                                      const Char_t * type    = "MC",
                                      const Char_t * period  = "LHC10H8",
                                      const Char_t * pass    = "PASS1",
                                      const Char_t * system  = "PBPB",
                                      const Double_t * initialParameters = 0x0,
                                      const Char_t * dedxtype= "",
                                      FitType = kSaturatedLund);
  static  TObjArray * GetResidualGraphs(THnSparseF * histPidQA, const Char_t * system, const Bool_t useV0s);
  static  TObjArray * GetResidualGraphsMC(THnSparseF * histPidQA, const Char_t * system);
  static  TObjArray * GetSeparation(THnSparseF * histPidQA, Int_t kParticle1, Int_t kParticle2);
  static  TObjArray * GetResponseFunctions(TF1* parametrisation, TObjArray* inputGraphs, const Char_t * type, const Char_t * period, const Char_t * pass, const Char_t * system, const Char_t * dedxtype);
  
  static TF1* SetUpFitFunction(const Double_t * initialParameters, AliTPCcalibResidualPID_Modified::FitType fitType, Float_t from, Float_t to, Bool_t isPPb, Bool_t isMC, Double_t* parametersBBForward);
  static void SetUpInputGraph(TGraphErrors* graphAll, Bool_t isMC, Bool_t useV0s);
  static TCanvas* CreateBBCanvas(TObjArray* inputGraphs, Bool_t isMC, TF1* func);
  static TCanvas* CreateResidualCanvas(TGraphErrors* graphAll, TF1* func);
  
  static  TF1*        FitBB(TObjArray* inputGraphs, Bool_t isMC, Bool_t isPPb, const Bool_t useV0s,
                            const Double_t * initialParameters = 0x0, FitType = kSaturatedLund);
  static Int_t MergeGraphErrors(TGraphErrors* mergedGraph, TCollection* li);
  
  static Double_t GetCutGeo() { return fgCutGeo; };
  static Double_t GetCutNcr() { return fgCutNcr; };
  static Double_t GetCutNcl() { return fgCutNcl; };
  
  static void SetCutGeo(Double_t value) { fgCutGeo = value; };
  static void SetCutNcr(Double_t value) { fgCutNcr = value; };
  static void SetCutNcl(Double_t value) { fgCutNcl = value; };
  
  static Bool_t TPCCutMIGeo(const AliVTrack* track, const AliVEvent* evt, TTreeStream* streamer = 0x0);
  static Bool_t TPCCutMIGeo(const AliVTrack* track, const AliInputEventHandler* evtHandler, TTreeStream* streamer = 0x0)
    { if (!evtHandler) return kFALSE; return TPCCutMIGeo(track, evtHandler->GetEvent(), streamer); };

  static TString GetStringFitType(Int_t fitType);
    
  protected:
  static Double_t fgCutGeo;  // Cut variable for TPCCutMIGeo concerning geometry
  static Double_t fgCutNcr;  // Cut variable for TPCCutMIGeo concerning num crossed rows
  static Double_t fgCutNcl;  // Cut variable for TPCCutMIGeo concerning num clusters
  
  static Double_t Lund(Double_t* xx, Double_t* par);
  static Double_t SaturatedLund(Double_t* xx, Double_t* par);
  static Double_t Aleph(Double_t* xx, Double_t* par);
  
  static void BinLogAxis(THnSparseF *h, Int_t axisNumber);
  static THnSparseF* InitialisePIDQAHist(TString name, TString title, Bool_t IsPbPb = kFALSE);
  static void  SetAxisNamesFromTitle(const THnSparseF *h);

  static void FitSlicesY(TH2 *hist, Double_t heightFractionForRange, Int_t cutThreshold, TString fitOption, TObjArray *arr);

  void FillV0PIDlist(AliESDEvent* esdEvent = 0x0);
  void ClearV0PIDlist();
  
  private:
  //
  //
  AliESDEvent *fESD;                   //! ESD object
  AliMCEvent  *fMC;                    //! MC object
  TObjArray * fOutputContainer;        //! output data container
  AliESDtrackCuts * fESDtrackCuts;     // basic cut variables for all non-V0 tracks
  AliESDtrackCuts * fESDtrackCutsV0;   // basic cut variables for all V0 tracks
  AliPIDResponse* fPIDResponse;        //! PID handling
  //
  
  Short_t fNumEtaCorrReqErrorsIssued;  // Number of times the error about eta correction issues have been displayed
  Short_t fNumMultCorrReqErrorsIssued; // Number of times the error about multiplicity correction issues have been displayed
  
  Bool_t fUseTPCCutMIGeo;   // Use geometrical cut for TPC
  
  Bool_t fUseMCinfo;         // Use MC info, if available
  
  Bool_t fIsPbpOrpPb;      // Pbp/pPb collision or something else?
  Bool_t fIsPbPb;          // PbPb collision?
  Double_t fZvtxCutEvent;  // Vertex z cut for the event (cm)
  
  AliESDv0KineCuts *fV0KineCuts;       //! ESD V0 kine cuts
  AliAnalysisUtils *fAnaUtils; //! Object to use analysis utils like pile-up rejection
  Bool_t fCutOnProdRadiusForV0el;      // Cut on production radius for V0 electrons
  Int_t fNumTagsStored;     // Number of entries of fV0tags
  Char_t* fV0tags;         //! Pointer to array with tags for identified particles from V0 decays
  Int_t* fV0motherIndex;   //! Pointer to array with index of the mother V0
  Int_t* fV0motherPDG;     //! Pointer to array with pdg of the mother V0

  Bool_t fProduceAllPadTypes, fProduceGlobal, fProduceShortPads, fProduceMediumPads, fProduceLongPads,fProduceOroc;
  THnSparseF * fHistPidQA;             //! histogram for the QA of the PID
  THnSparseF * fHistPidQAshort;        //! histogram for the QA of the PID short pads
  THnSparseF * fHistPidQAmedium;       //! histogram for the QA of the PID med pads
  THnSparseF * fHistPidQAlong;         //! histogram for the QA of the PID long pads
  THnSparseF * fHistPidQAoroc;         //! histogram for the QA of the PID full oroc
  //
  Bool_t fProduceTPCSignalSparse;            // for setter
  Bool_t fCorrectdEdxEtaDependence;          // Correct eta dependence for fHistPidQA (NOTE: Not done for the pad-specific THnSparses)
  Bool_t fCorrectdEdxMultiplicityDependence; // Correct multiplicity dependence for fHistPidQA (NOTE: Not done for the pad-specific THnSparses)
  Bool_t fCorrectdEdxPileupDependence;       // Correct pileup dependence for fHistPidQA (NOTE: Not done for the pad-specific THnSparses)
  //
  THnSparseF * fThnspTpc;              //! thnsparse containing the data
  //
  //
  
  Bool_t fWriteAdditionalOutput; // Also fill histos/trees for QA etc. and write them
  
  // QA histos
  TObjArray* fQAList;           //! Array with QA histos
  TH1F* fhInvMassGamma;         //! Histogram with inv. mass of gamma
  TH1F* fhInvMassK0s;           //! Histogram with inv. mass of K0s
  TH1F* fhInvMassLambda;        //! Histogram with inv. mass of lambda
  TH1F* fhInvMassAntiLambda;    //! Histogram with inv. mass of anti-lambda
  
  TH2F* fhArmenterosAll;        //! Histogram with armenteros plot for all V0s
  TH2F* fhArmenterosGamma;      //! Histogram with armenteros plot for gamma
  TH2F* fhArmenterosK0s;        //! Histogram with armenteros plot for K0s
  TH2F* fhArmenterosLambda;     //! Histogram with armenteros plot for lambda
  TH2F* fhArmenterosAntiLambda; //! Histogram with armenteros plot for anti-lambda
  
  // QA histos for shared clusters
  THnSparseF* fHistSharedClusQAV0Pi;  //! Histogram with shared clusters QA for V0 pi
  THnSparseF* fHistSharedClusQAV0Pr;  //! Histogram with shared clusters QA for V0 pr
  THnSparseF* fHistSharedClusQAV0El;  //! Histogram with shared clusters QA for V0 el
    
  
  // TTree stuff for advanced studies (local track density, ...)
  TTree* fTreeV0El;               //! Tree with V0 el and closest neighbour tracks and V0 sisters
  TTree* fTreeV0Pi;               //! Tree with V0 pi and closest neighbour tracks and V0 sisters
  TTree* fTreeV0Pr;               //! Tree with V0 pr and closest neighbour tracks and V0 sisters
  Double_t fTree_dEdx_tr;         //! Tree: dEdx of track
  Double_t fTree_dEdx_nb;         //! Tree: dEdx of neighbour
  Double_t fTree_dEdx_vs;         //! Tree: dEdx of V0 sister
  Double_t fTree_dEdxExpected_tr; //! Tree: dEdx_expected of track
  Double_t fTree_p_TPC_tr;        //! Tree: TPC momentum of track
  Double_t fTree_p_TPC_nb;        //! Tree: TPC momentum of neighbour
  Double_t fTree_p_TPC_vs;        //! Tree: TPC momentum of V0 sister
  Double_t fTree_BtimesChargeOverPt_tr; //! Tree: mag field times charge/pT of track
  Double_t fTree_BtimesChargeOverPt_nb; //! Tree: mag field times charge/pT of neighbour
  Double_t fTree_BtimesChargeOverPt_vs; //! Tree: mag field times charge/pT of V0 sister
  Double_t fTree_tanTheta_tr;     //! Tree: tan(theta) of track
  Double_t fTree_tanTheta_nb;     //! Tree: tan(theta) of neighbour
  Double_t fTree_tanTheta_vs;     //! Tree: tan(theta) of V0 sister
  Double_t fTree_distance_nb;     //! Tree: distance on TPC cylinder of track and neighbour
  Double_t fTree_distance_vs;     //! Tree: distance on TPC cylinder of track and V0 sister
    
    
    
    
    
    
    
    //Reduced Tree (Alberto Caliva)
    TTree *fReducedTree;//!
    
    //Tree Variables (Alberto Caliva)
    Int_t    fParticleIDcode;//!
    Int_t    fNumberOfTPCclusters_dEdx;//!
    Double_t fTPC_dEdx_au;//!
    Double_t fPseudorapidity;//!
    Double_t fCentralityPercentile;//!
    Double_t fMomentumTPCinnnerWall;//!
    Double_t fNumberOfSigmasTOF_Elec;//!
    Double_t fNumberOfSigmasTOF_Pion;//!
    Double_t fNumberOfSigmasTOF_Kaon;//!
    Double_t fNumberOfSigmasTOF_Prot;//!
    

    
  
  AliTPCcalibResidualPID_Modified(const AliTPCcalibResidualPID_Modified&); // not implemented
  AliTPCcalibResidualPID_Modified& operator=(const AliTPCcalibResidualPID_Modified&); // not implemented
  
  ClassDef(AliTPCcalibResidualPID_Modified, 6);
};

#endif
