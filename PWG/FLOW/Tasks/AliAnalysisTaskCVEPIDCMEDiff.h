#ifndef AliAnalysisTaskCVEPIDCMEDiff_cxx
#define AliAnalysisTaskCVEPIDCMEDiff_cxx
#include <vector>
#include <map>
#include <unordered_map>
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"


class AliAnalysisTaskCVEPIDCMEDiff : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskCVEPIDCMEDiff();
  AliAnalysisTaskCVEPIDCMEDiff(const char* name);
  virtual ~AliAnalysisTaskCVEPIDCMEDiff();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  // Switch
  void IfCalculateLambdaHadron(bool bCalculateLambdaHadron) { this->isCalculateLambdaHadron = bCalculateLambdaHadron; }
  void IfDebug(bool bDebug) { this->fDebug = bDebug; }
  void IfDoNUE(bool bDoNUE) { this->isDoNUE = bDoNUE; }
  void IfDoLambdaNUE(bool bDoLambdaNUE) { this->isDoLambdaNUE = bDoLambdaNUE; }
  void IfV0DaughterUseTOF(bool bV0DaughterUseTOF) { this->isV0DaughterUseTOF = bV0DaughterUseTOF; }
  void IfNarrowDcaCuts768(bool bNarrowDcaCuts768) { this->isNarrowDcaCuts768 = bNarrowDcaCuts768; }
  void IfProtonCustomizedDCACut(bool bProtonCustomizedDCACut) { this->isProtonCustomizedDCACut = bProtonCustomizedDCACut; }
  void IfUsePionRejection(bool bUsePionRejection) { this->isUsePionRejection = bUsePionRejection; }
  void IfTightPileUp(bool bTightPileUp) { this->isTightPileUp = bTightPileUp; }

  // read in
  void SetListForNUE(TList* flist) { this->fListNUE = (TList*)flist->Clone(); }
  void SetListForNUA(TList* flist) { this->fListNUA = (TList*)flist->Clone(); }
  void SetListForVZEROCalib(TList* flist) { this->fListVZEROCalib = (TList*)flist->Clone(); }

  // Global
  void SetTrigger(TString trigger) { this->fTrigger = trigger; }
  void SetPeriod(TString period) { this->fPeriod = period; }
  // Event
  void SetVzCut(double vzCut) { this->fVzCut = vzCut; }
  void SetCentCut(float centDiffCut) { this->fCentDiffCut = centDiffCut; }
  // Plane
  void SetPlaneEstimator(TString planeEstimator) { this->fPlaneEstimator = planeEstimator; }
  void SetPlanePtMin(float planePtMin) { this->fPlanePtMin = planePtMin; }
  void SetPlanePtMax(float planePtMax) { this->fPlanePtMax = planePtMax; }
  // Track
  void SetFilterBit(int filterBit) { this->fFilterBit = filterBit; }
  void SetNclsCut(int nclsCut) { this->fNclsCut = nclsCut; }
  void SetChi2Max(float chi2Max) { this->fChi2Max = chi2Max; }
  void SetChi2Min(float chi2Min) { this->fChi2Min = chi2Min; }
  void SetDCAcutXY(float dcaCutxy) { this->fDcaCutXY = dcaCutxy; }
  void SetDCAcutZ(float dcaCutz) { this->fDcaCutZ = dcaCutz; }
  void SetPtMin(float ptMin) { this->fPtMin = ptMin; }
  void SetPtMax(float ptMax) { this->fPtMax = ptMax; }
  void SetEtaCut(float etaCut) { this->fEtaCut = etaCut; }
  void SetDedxCut(float dedxCut) { this->fDedxCut = dedxCut; }
  void SetProtonPtMin(double protonPtMin) { this->fProtonPtMin = protonPtMin; }
  void SetProtonPtMax(double protonPtMax) { this->fProtonPtMax = protonPtMax; }
  void SetAntiProtonPtMin(double antiprotonPtMin) { this->fAntiProtonPtMin = antiprotonPtMin; }
  void SetAntiProtonPtMax(double antiprotonPtMax) { this->fAntiProtonPtMax = antiprotonPtMax; }
  // PID
  void SetNSigmaTPCCut(double nSigmaTPC) { this->fNSigmaTPCCut = nSigmaTPC; }
  void SetNSigmaTOFCut(double nSigmaTOF) { this->fNSigmaTOFCut = nSigmaTOF; }
  // V0
  void SetV0CPAMin(double v0CPAMin) { this->fV0CPAMin = v0CPAMin; }
  void SetV0DecayLengthMax(double v0DecayLengthMax) { this->fV0DecayLengthMax = v0DecayLengthMax; }
  void SetV0DecayLengthMin(double v0DecayLengthMin) { this->fV0DecayLengthMin = v0DecayLengthMin; }
  void SetV0DCAToPrimVtxMax(double v0DCAToPrimVtxMax) { this->fV0DCAToPrimVtxMax = v0DCAToPrimVtxMax; }
  void SetV0DcaBetweenDaughtersMax(double v0DcaBetweenDaughtersMax) { this->fV0DcaBetweenDaughtersMax = v0DcaBetweenDaughtersMax; }
  // V0 Daughter Cut
  void SetDaughtersPtMax(double daughtersPtMax) { this->fDaughtersPtMax = daughtersPtMax; }
  void SetDaughtersEtaMax(double daughtersEtaMax) { this->fDaughtersEtaMax = daughtersEtaMax; }
  void SetDaughtersTPCNclsMin(double daughtersTPCNclsMin) { this->fDaughtersTPCNclsMin = daughtersTPCNclsMin; }
  void SetDaughtersDCAToPrimVtxMin(double daughtersDCAToPrimVtxMin) { this->fDaughtersDCAToPrimVtxMin = daughtersDCAToPrimVtxMin; }
  void SetRatioCrossedRowsFindable(double ratioCrossedRowsFindable) { this->fRatioCrossedRowsFindable = ratioCrossedRowsFindable; }
  void SetPosProtonTPCNsigma(float V0PosProtonTPCNsigma) { this->fV0PosProtonTPCNsigma = V0PosProtonTPCNsigma; }
  void SetNegPionTPCNsigma(float V0NegPionTPCNsigma) { this->fV0NegPionTPCNsigma = V0NegPionTPCNsigma; }
  void SetNegProtonTPCNsigma(float V0NegProtonTPCNsigma) { this->fV0NegProtonTPCNsigma = V0NegProtonTPCNsigma; }
  void SetPosPionTPCNsigma(float V0PosPionTPCNsigma) { this->fV0PosPionTPCNsigma = V0PosPionTPCNsigma; }
  void SetPosProtonTOFNsigma(float V0PosProtonTOFNsigma) { this->fV0PosProtonTOFNsigma = V0PosProtonTOFNsigma; }
  void SetNegPionTOFNsigma(float V0NegPionTOFNsigma) { this->fV0NegPionTOFNsigma = V0NegPionTOFNsigma; }
  void SetNegProtonTOFNsigma(float V0NegProtonTOFNsigma) { this->fV0NegProtonTOFNsigma = V0NegProtonTOFNsigma; }
  void SetPosPionTOFNsigma(float V0PosPionTOFNsigma) { this->fV0PosPionTOFNsigma = V0PosPionTOFNsigma; }
  // Lambda Cut
  void SetLambdaPtMin(double lambdaPtMin) { this->fLambdaPtMin = lambdaPtMin; }
  void SetLambdaPtMax(double lambdaPtMax) { this->fLambdaPtMax = lambdaPtMax; }
  void SetAntiLambdaPtMin(double antiLambdaPtMin) { this->fAntiLambdaPtMin = antiLambdaPtMin; }
  void SetAntiLambdaPtMax(double antiLambdaPtMax) { this->fAntiLambdaPtMax = antiLambdaPtMax; }
  void SetLambdaMassLeftCut(double lambdaMassRightCut) { this->fLambdaMassRightCut = lambdaMassRightCut; }
  void SetLambdaMassRightCut(double lambdaMassLeftCut) { this->fLambdaMassLeftCut = lambdaMassLeftCut; }
  void SetAntiLambdaMassLeftCut(double antiLambdaMassRightCut) { this->fAntiLambdaMassRightCut = antiLambdaMassRightCut; }
  void SetAntiLambdaMassRightCut(double antiLambdaMassLeftCut) { this->fAntiLambdaMassLeftCut = antiLambdaMassLeftCut; }
  void SetNMassBins(int nMassBins) { this->fNMassBins = nMassBins; }

 private:
  ////////////////////////
  // Procedural function
  ////////////////////////
  double GetTPCPlane();
  double GetV0CPlane();
  void ResetVectors();
  bool LoopTracks();
  bool LoopV0s();
  bool PairV0Trk();

  ////////////////////////
  // Functional function
  ////////////////////////
  // Read in
  bool LoadCalibHistForThisRun(); //deal with all the readin
  // Pile-up
  bool RejectEvtTFFit();
  // Track
  bool AcceptAODTrack(AliAODTrack* track);
  bool CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck);
  double GetNUECor(int charge, double pt);
  double GetPIDNUECor(int pdgcode, double pt);
  double GetNUACor(int charge, double phi, double eta, double vz);
  // V0
  bool IsGoodV0(AliAODv0* aodV0);
  bool IsGoodDaughterTrack(const AliAODTrack* track);
  int GetLambdaCode(const AliAODTrack* pTrack, const AliAODTrack* ntrack);
  // Plane
  double GetTPCPlaneNoAutoCorr(std::vector<int> vec_id);
  inline double GetEventPlane(double qx, double qy, double harmonic);
  // Get DCA
  bool GetDCA(double &dcaxy, double &dcaz, AliAODTrack* track);
  // Sum pT bin, Delta eta bin
  inline double GetSumPtBin(double sumPt);
  inline double GetDeltaEtaBin(double deltaEta);
  inline double GetDCABin(double dca);

  //////////////////////
  // Switch           //
  //////////////////////
  bool isCalculateLambdaHadron;
  bool isDoNUE;
  bool isDoLambdaNUE;
  bool isV0DaughterUseTOF;
  bool isNarrowDcaCuts768;
  bool isProtonCustomizedDCACut;
  bool isUsePionRejection;
  bool isTightPileUp;

  //////////////////////
  // Cuts and options //
  //////////////////////
  // Global
  TString fTrigger; //
  TString fPeriod;  // period
  // Event
  double fVzCut;      // vz cut
  float fCentDiffCut; // centrality restriction for V0M and TRK
  // Plane
  TString fPlaneEstimator; // TPC or V0C
  float fPlanePtMin;
  float fPlanePtMax;
  // Track
  int fFilterBit;          // AOD filter bit selection
  int fNclsCut;            // ncls cut for all tracks
  float fChi2Max;          // upper limmit for chi2
  float fChi2Min;          // lower limmit for chi2
  float fDcaCutXY;            // dcaxy cut for all tracks
  float fDcaCutZ;             // dcaz cut for all tracks
  float fPtMin;            // minimum pt for tracks
  float fPtMax;            // maximum pt for tracks
  float fEtaCut;           // eta cut for tracks
  float fDedxCut;          // dedx cut for tracks
  double fProtonPtMin;     // Min pt for proton
  double fProtonPtMax;     // Max pt for proton
  double fAntiProtonPtMin; // Min pt for anti-proton
  double fAntiProtonPtMax; // Max pt for anti-proton
  // PID
  float fNSigmaTPCCut;
  float fNSigmaTOFCut;
  // V0
  double fV0CPAMin;                 //
  double fV0DecayLengthMin;         //
  double fV0DecayLengthMax;         //
  double fV0DCAToPrimVtxMax;        //
  double fV0DcaBetweenDaughtersMax; //
  // V0 Daughter
  double fDaughtersPtMax;           //
  double fDaughtersEtaMax;          //
  double fDaughtersTPCNclsMin;      //
  double fDaughtersDCAToPrimVtxMin; //
  double fRatioCrossedRowsFindable; //
  float fV0PosProtonTPCNsigma;      //
  float fV0NegPionTPCNsigma;        //
  float fV0NegProtonTPCNsigma;      //
  float fV0PosPionTPCNsigma;        //
  float fV0PosProtonTOFNsigma;      //
  float fV0NegPionTOFNsigma;        //
  float fV0NegProtonTOFNsigma;      //
  float fV0PosPionTOFNsigma;        //
  // Lambda Cut
  double fLambdaPtMin;              //
  double fLambdaPtMax;              //
  double fAntiLambdaPtMin;          //
  double fAntiLambdaPtMax;          //
  double fLambdaRapidityMax;        //
  double fAntiLambdaRapidityMax;    //
  double fLambdaMassRightCut;       //
  double fLambdaMassLeftCut;        //
  double fAntiLambdaMassRightCut;   //
  double fAntiLambdaMassLeftCut;    //
  double fLambdaMassMean;           //
  // Number of Mass Bins
  int fNMassBins;

  ///////////////////The following files are from the data//////////////////////////////////
  /////////////
  // Handles //
  /////////////
  AliAODEvent* fAOD;            // aod Event
  AliPIDResponse* fPIDResponse; // PID Handler
  AliAnalysisUtils* fUtils;     // Event Selection Options
  AliMultSelection* fMultSel;

  ////////////////////////////////
  // Global Variables from data //
  ////////////////////////////////
  double fVertex[3]; // vetex
  int fRunNum;       // runnumber
  int fOldRunNum;    // latest runnumber
  int fRunNumBin;    // runnumer bin; 10:139510...; 11:170387...; 15HIR:246994...
  int fVzBin;        // vertex z bin
  float fCent;       // centrality
  int fCentBin;      // centrality bin: 0-10
  float fCentV0M;    // Centrality V0M
  float fCentTRK;    // Centrality TRK
  float fCentSPD0;   // Centrality SPD0
  float fCentSPD1;   // Centrality SPD1
  // Variable to get TPC Plane
  double fSumQ2xTPC;
  double fSumQ2yTPC;
  double fWgtMultTPC;
  // Plane
  double fPsi2;

  // Plane tracks Map key:id value:(phi,weight)
  std::unordered_map<int, std::vector<double>> mapTPCTrksIDPhiWgt;
  
  // Vector for particles from Tracks [pt,eta,phi,id,pdgcode,weight,pidweight,dcaxy]
  std::vector<std::array<double,8>> vecParticle;
  // Vector for V0s [pt,eta,phi,id,pdgcode,weight,mass,id1,id2]
  std::vector<std::array<double,9>> vecParticleV0;

  ///////////////////The following files are read from external sources////////////////////
  ////////////////////////
  // Pile up Function
  ////////////////////////
  TF1* fSPDCutPU;
  TF1* fV0CutPU;
  TF1* fCenCutLowPU;
  TF1* fCenCutHighPU;
  TF1* fMultCutPU;
  ////////////////////////
  // NUE
  ////////////////////////
  TList* fListNUE; // read list for NUE
  TH1D* hNUEweightPlus;
  TH1D* hNUEweightMinus;
  ////////////////////////
  // NUA
  ////////////////////////
  TList* fListNUA; // read lists for NUA
  TH3F* hCorrectNUAPos;
  TH3F* hCorrectNUANeg;
  ////////////////////////
  // VZERO
  ////////////////////////
  TList* fListVZEROCalib; // read list for V0 Calib
  //for recenter
  AliOADBContainer* contQxncm;
  AliOADBContainer* contQyncm;
  TH1D* hQx2mV0C;
  TH1D* hQy2mV0C;
  // for gain equalization
  TH2F* fHCorrectV0ChWeghts;


  ///////////////////The following files will be saved//////////////////////////////////
  //////////////
  // QA Plots //
  //////////////
  TList* fQAList;
  // General QA
  // Event-wise
  TH1D* fEvtCount;
  std::map<int, int>* runNumList;
  TH1I* fHistRunNumBin;
  TH1D* fHistCent[2];
  TH1D* fHistVz[2];
  TH2D* fHist2CentQA[8];
  TH2D* fHist2MultCentQA[2]; // need for cheak pile up
  TH2D* fHist2MultMultQA[6]; // need for cheak pile up
  // Track-wise
  TH1D* fHistPt;
  TH1D* fHistEta;
  TH1D* fHistNhits;
  TH2D* fHist2PDedx;
  TH1D* fHistDcaXY;
  TH1D* fHistDcaZ;
  TH1D* fHistPhi[2];
  

  //Proton QA
  TH1D* fHistProtonPt;
  TH1D* fHistProtonEta;
  TH1D* fHistProtonPhi;
  TH1D* fHistProtonDcaXY;
  TH1D* fHistProtonDcaZ;
  TH2D* fHist2ProtonSigTPC;
  TH2D* fHist2ProtonSigTOF;
  TH1D* fHistAntiProtonPt;
  TH1D* fHistAntiProtonEta;
  TH1D* fHistAntiProtonPhi;
  TH1D* fHistAntiProtonDcaXY;
  TH1D* fHistAntiProtonDcaZ;
  TH2D* fHist2AntiProtonSigTPC;
  TH2D* fHist2AntiProtonSigTOF;

  // V0s QA
  TH1D* fHistV0Pt;              // !Raw V0s' pT
  TH1D* fHistV0Eta;             // !Raw V0s' eta
  TH1D* fHistV0DcatoPrimVertex; // !Raw V0s' DcatoPV
  TH1D* fHistV0CPA;             // !Raw V0s' CPA(cosine pointing angle)
  TH1D* fHistV0DecayLength;     // !Raw V0s' DecayLength
  TH1D* fHistV0NegDaughterDca;  // !Raw V0s' NegDaughterDca
  TH1D* fHistV0PosDaughterDca;  // !Raw V0s' PosDaughterDca
  // Lambda QA
  //[0]:Before the Mass Cut [1]:After the Mass Cut
  TH1D* fHistLambdaPt[2];                  //
  TH1D* fHistLambdaEta[2];                 //
  TH1D* fHistLambdaPhi[2];                 //
  TH1D* fHistLambdaDcaToPrimVertex[2];     //
  TH1D* fHistLambdaNegDaugtherDca[2];      //
  TH1D* fHistLambdaPosDaugtherDca[2];      //
  TH1D* fHistLambdaCPA[2];                 //
  TH1D* fHistLambdaDecayLength[2];         //
  TH3D* fHist3LambdaCentPtMass[2];         //
  TH2D* fHist2LambdaMassPtY[2];            //
  TH1D* fHistAntiLambdaPt[2];              //
  TH1D* fHistAntiLambdaEta[2];             //
  TH1D* fHistAntiLambdaPhi[2];             //
  TH1D* fHistAntiLambdaDcaToPrimVertex[2]; //
  TH1D* fHistAntiLambdaNegDaugtherDca[2];  //
  TH1D* fHistAntiLambdaPosDaugtherDca[2];  //
  TH1D* fHistAntiLambdaCPA[2];             //
  TH1D* fHistAntiLambdaDecayLength[2];     //
  TH3D* fHist3AntiLambdaCentPtMass[2];     //
  TH2D* fHist2AntiLambdaMassPtY[2];        //


  /////////////
  // Results //
  /////////////
  TList* fResultsList;
  // Plane
  TH2D* fHist2Psi2;

  // Inv Mass 
  TH2D* fHist2LambdaMass[4];           //![0]:Λ    [1]:Λ      [2]:Λbar    [3]:Λbar
  TH3D* fHist3LambdaProtonMassSPt[4];  //![0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  TH3D* fHist3LambdaProtonMassDEta[4]; //![0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar

  // Diff δ(ΔpT)
  TProfile3D* fProfile3DDiffDeltaLambdaProtonMassSPt[4]; //![0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  // Diff δ(Δη)
  TProfile3D* fProfile3DDiffDeltaLambdaProtonMassDEta[4]; //![0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  // Diff γ(SpT)
  TProfile3D* fProfile3DDiffGammaLambdaProtonMassSPt[4]; //![0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  // Diff γ(Δη)
  TProfile3D* fProfile3DDiffGammaLambdaProtonMassDEta[4]; //![0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar

  TProfile2D* fProfile2DDeltaLambdaHadronMass[4]; //![0]:Λ-h+ [1]:Λ-h-   [2]:Λbar-h+ [3]:Λbar-h-
  TProfile2D* fProfile2DGammaLambdaHadronMass[4]; //![0]:Λ-h+ [1]:Λ-h-   [2]:Λbar-h+ [3]:Λbar-h-


  TH3D* fHist3LambdaProtonMassDCA[4];
  //Diff δ(dca)
  TProfile3D* fProfile3DDiffDeltaLambdaProtonMassDCA[4]; //![0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar
  //Diff γ(dca)
  TProfile3D* fProfile3DDiffGammaLambdaProtonMassDCA[4]; //![0]:Λ-p  [1]:Λ-pbar [2]:Λbar-p  [3]:Λbar-pbar

  AliAnalysisTaskCVEPIDCMEDiff(const AliAnalysisTaskCVEPIDCMEDiff&);
  AliAnalysisTaskCVEPIDCMEDiff& operator=(const AliAnalysisTaskCVEPIDCMEDiff&);

  ClassDef(AliAnalysisTaskCVEPIDCMEDiff, 1);
};

#endif

