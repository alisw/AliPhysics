#ifndef AliAnalysisTaskPIDCME_cxx
#define AliAnalysisTaskPIDCME_cxx
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

class AliAnalysisTaskPIDCME : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskPIDCME();
  AliAnalysisTaskPIDCME(const char* name);
  virtual ~AliAnalysisTaskPIDCME();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  // Switch for adjustion in AddTsk.C in train parameter setting
  void IfDebug(bool bDebug) { this->fDebug = bDebug; }
  void IfUseTPCPlane(bool bUseTPCPlane) { this->isUseTPCPlane = bUseTPCPlane; }
  void IfUseVZEROPlane(bool bUseVZEROPlane) { this->isUseVZEROPlane = bUseVZEROPlane; }
  void IfUseZDCPlane(bool bUseZDCPlane) { this->isUseZDCPlane = bUseZDCPlane; }
  void IfDoNUE(bool bDoNUE) { this->isDoNUE = bDoNUE; }
  void IfDoNUA(bool bDoNUA) { this->isDoNUA = bDoNUA; }
  void IfQATPC(bool bQATPC) { this->isQATPC = bQATPC; }
  void IfQAVZERO(bool bQAVZERO) { this->isQAVZERO = bQAVZERO; }
  void IfQAZDC(bool bQAZDC) { this->isQAZDC = bQAZDC; }
  void IfNarrowDcaCuts768(bool bNarrowDcaCuts768) { this->isNarrowDcaCuts768 = bNarrowDcaCuts768; }
  void IfProtonCustomizedDCACut(bool bProtonCustomizedDCACut) { this->isProtonCustomizedDCACut = bProtonCustomizedDCACut; }
  void IfUsePionRejection(bool bUsePionRejection) { this->isUsePionRejection = bUsePionRejection; }

  void IfTightPileUp(bool bTightPileUp) { this->isTightPileUp = bTightPileUp; }

  void IfCalculatePIDFlow(bool bCalculatePIDFlow) { this->isCalculatePIDFlow = bCalculatePIDFlow; }
  void IfCalculateDiffResult(bool bCalculateDiffResult) { this->isCalculateDiffResult = bCalculateDiffResult; }
  void IfCalculateDeltaPhiSumPhi(bool bCalculateDeltaPhiSumPhi) { this->isCalculateDeltaPhiSumPhi = bCalculateDeltaPhiSumPhi; }

  void IfCalculatePionKaon(bool bCalculatePionKaon) {this->isCalculatePionKaon = bCalculatePionKaon; }
  void IfCalculatePionProton(bool bCalculatePionProton) {this->isCalculatePionProton = bCalculatePionProton; }
  void IfCalculateKaonProton(bool bCalculateKaonProton) {this->isCalculateKaonProton = bCalculateKaonProton; }
  void IfCalculatePionPion(bool bCalculatePionPion) {this->isCalculatePionPion = bCalculatePionPion; }
  void IfCalculateKaonKaon(bool bCalculateKaonKaon) {this->isCalculateKaonKaon = bCalculateKaonKaon; }
  void IfCalculateProtonProton(bool bCalculateProtonProton) {this->isCalculateProtonProton = bCalculateProtonProton; }
  void IfCalculateHadronHadron(bool bCalculateHadronHadron) {this->isCalculateHadronHadron = bCalculateHadronHadron; }
  
  void IfUseOneSideTPCPlane(bool bUseOneSideTPCPlane) { this->isUseOneSideTPCPlane = bUseOneSideTPCPlane; }
  void IfOpenPIDSingletrk(bool bOpenPIDSingletrk) { this->isOpenPIDSingletrk = bOpenPIDSingletrk; }
  void IfOpenSsandOsSelfCheck(bool bOpenSsandOsSelfCheck) {this->isOpenSsandOsSelfCheck = bOpenSsandOsSelfCheck; }

  // read in
  void SetListForNUE(TList* flist) { this->fListNUE = (TList*)flist->Clone(); }
  void SetListForNUA(TList* flist) { this->fListNUA = (TList*)flist->Clone(); }
  void SetListForVZEROCalib(TList* flist) { this->fListVZEROCalib = (TList*)flist->Clone(); }
  void SetListForZDCCalib(TList* flist) { this->fListZDCCalib = (TList*)flist->Clone(); }

  // Global
  void SetTrigger(TString trigger) { this->fTrigger = trigger; }
  void SetPeriod(TString period) { this->fPeriod = period; }
  // Event
  void SetVzCut(double vzCut) { this->fVzCut = vzCut; }
  void SetCentCut(float centDiffCut) { this->fCentDiffCut = centDiffCut; }
  // Plane
  void SetPlanePtMin(float planePtMin) { this->fPlanePtMin = planePtMin; }
  void SetPlanePtMax(float planePtMax) { this->fPlanePtMax = planePtMax; }
  void SetPlaneEtaGapPos(float etaGapPos) { this->fEtaGapPos = etaGapPos; }
  void SetPlaneEtaGapNeg(float etaGapNeg) { this->fEtaGapNeg = etaGapNeg; }
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
  void SetPionPtMin(double PionPtMin) { this->fPionPtMin = PionPtMin; }
  void SetPionPtMax(double PionPtMax) { this->fPionPtMax = PionPtMax; }
  void SetKaonPtMin(double KaonPtMin) { this->fKaonPtMin = KaonPtMin; }
  void SetKaonPtMax(double KaonPtMax) { this->fKaonPtMax = KaonPtMax; }
  void SetAntiProtonPtMin(double antiprotonPtMin) { this->fAntiProtonPtMin = antiprotonPtMin; }
  void SetAntiProtonPtMax(double antiprotonPtMax) { this->fAntiProtonPtMax = antiprotonPtMax; }
  void SetAntiPionPtMin(double antiPionPtMin) { this->fAntiPionPtMin = antiPionPtMin; }
  void SetAntiPionPtMax(double antiPionPtMax) { this->fAntiPionPtMax = antiPionPtMax; }
  void SetAntiKaonPtMin(double antiKaonPtMin) { this->fAntiKaonPtMin = antiKaonPtMin; }
  void SetAntiKaonPtMax(double antiKaonPtMax) { this->fAntiKaonPtMax = antiKaonPtMax; }
  // PID
  void SetPtMinforRMSPion(double PtMinforRMSPion) { this->fPtMinforRMSPion = PtMinforRMSPion; }
  void SetPtMinforRMSKaon(double PtMinforRMSKaon) { this->fPtMinforRMSKaon = PtMinforRMSKaon; }
  void SetPtMinforRMSProton(double PtMinforRMSProton) { this->fPtMinforRMSProton = PtMinforRMSProton; }
  void SetNSigmaTPCCutPion(double NSigmaTPCCutPion) { this->fNSigmaTPCCutPion = NSigmaTPCCutPion; }
  void SetNSigmaTPCCutKaon(double NSigmaTPCCutKaon) { this->fNSigmaTPCCutKaon = NSigmaTPCCutKaon; }
  void SetNSigmaTPCCutProton(double NSigmaTPCCutProton) { this->fNSigmaTPCCutProton = NSigmaTPCCutProton; }
  // ESE Qn
  void SetNQ2Bins(int nQ2Bins) { this->fNQ2Bins = nQ2Bins; }

 private:
  ////////////////////////
  // Procedural function
  ////////////////////////
  bool GetTPCPlane();
  bool GetVZEROPlane();
  bool GetZDCPlane();
  bool GetZDCPlaneLsFit();
  void ResetVectors();
  bool LoopTracks();
  bool PairTrkTrk();

  ////////////////////////
  // Functional function
  ////////////////////////
  // Read in
  bool LoadCalibHistForThisRun(); //deal with all the readin
  // Pile-up
  bool RejectEvtMultComp();
  bool RejectEvtTFFit();
  bool RejectEvtTPCITSfb32TOF();
  bool AODPileupCheck();
  bool PileUpMultiVertex();
  bool RemovalForRun1();
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  // Track
  bool AcceptAODTrack(AliAODTrack* track);
  bool CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck);
  double GetNUECor(int charge, double pt);
  double GetPIDNUECor(int pdgcode, double pt);
  double GetNUACor(int charge, double phi, double eta, double vz);
  int DetermineTheBestTPCPID(AliAODTrack* track, int PIDPoints);
  // Plane
  std::array<double,2> GetOneSideTPCPlaneNoAutoCorr(std::vector<int> vec_id);
  double GetTPCPlaneNoAutoCorr(std::vector<int> vec_id);
  inline double GetEventPlane(double qx, double qy, double harmonic);
  // Range phi
  inline double RangeDPhi(double dphi);
  // Get DCA
  bool GetDCA(double &dcaxy, double &dcaz, AliAODTrack* track);
  // Get q2 Bin
  double GetQ2Bin(double q2);

  //////////////////////
  // Switch           //
  //////////////////////
  bool isUseTPCPlane;
  bool isUseVZEROPlane;
  bool isUseZDCPlane;
  bool isDoNUE;
  bool isDoNUA;
  bool isQATPC;
  bool isQAVZERO; // flag for V0 qn QA
  bool isQAZDC;
  bool isNarrowDcaCuts768;
  bool isProtonCustomizedDCACut;
  bool isUsePionRejection;

  bool isCalculatePIDFlow;        // if fill PID Flow
  bool isCalculateDiffResult;     // if fill Diff Reslut
  bool isCalculateDeltaPhiSumPhi; // if fill C(delta_phi)

  bool isCalculatePionKaon;
  bool isCalculatePionProton;
  bool isCalculateKaonProton;
  bool isCalculatePionPion;
  bool isCalculateKaonKaon;
  bool isCalculateProtonProton;
  bool isCalculateHadronHadron;
  
  bool isUseOneSideTPCPlane;

  bool isTightPileUp;
  bool isOpenPIDSingletrk;
  bool isOpenSsandOsSelfCheck;

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
  float fPlanePtMin;
  float fPlanePtMax;
  float fEtaGapPos; // value for the Eta Gap Pos
  float fEtaGapNeg;
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
  double fPionPtMin;    
  double fPionPtMax;
  double fKaonPtMin;
  double fKaonPtMax;
  double fProtonPtMin;     // Min pt for proton
  double fProtonPtMax;     // Max pt for proton

  double fAntiPionPtMin;    
  double fAntiPionPtMax;
  double fAntiKaonPtMin;
  double fAntiKaonPtMax;
  double fAntiProtonPtMin; // Min pt for anti-proton
  double fAntiProtonPtMax; // Max pt for anti-proton
  // PID
  float fPtMinforRMSPion;
  float fPtMinforRMSKaon;
  float fPtMinforRMSProton;
  float fNSigmaTPCCutPion;
  float fNSigmaRMSCutPion;
  float fNSigmaTPCCutKaon;
  float fNSigmaRMSCutKaon;
  float fNSigmaTPCCutProton;
  float fNSigmaRMSCutProton;
  // Number of Qn Bins for ESE 
  int fNQ2Bins;                     //
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
  double fSumQ2xTPCPos;
  double fSumQ2yTPCPos;
  double fWgtMultTPCPos;
  double fSumQ2xTPCNeg;
  double fSumQ2yTPCNeg;
  double fWgtMultTPCNeg;
  double fSumQ2xTPC;
  double fSumQ2yTPC;
  double fWgtMultTPC;
  // Plane
  double fPsi2TPCPos;
  double fPsi2TPCNeg;
  double fPsi2V0C;
  double fPsi2V0A;
  double fPsi1ZNC;
  double fPsi1ZNA;
  double fPsi2TPC;
  double fQ2Bin;
  // Do we get the right plane?
  bool isRightTPCPlane;
  bool isRightVZEROPlane;
  bool isRightZDCPlane;

  // Plane tracks Map key:id value:(phi,weight)
  std::unordered_map<int, std::vector<double>> mapTPCPosTrksIDPhiWgt;
  std::unordered_map<int, std::vector<double>> mapTPCNegTrksIDPhiWgt;
  std::unordered_map<int, std::vector<double>> mapTPCTrksIDPhiWgt;
  
  // Vector for particles from Tracks [pt,eta,phi,id,pdgcode,weight,pidweight]
  std::vector<std::array<double,7>> vecParticle;

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
  // 10h/15o
  TList* fListNUE; // read list for NUE
  TH1D* hNUEweightPlus;
  TH1D* hNUEweightMinus;
  ////////////////////////
  // NUA
  ////////////////////////
  TList* fListNUA; // read lists for NUA
  // 10h
  TH2D* hNUAweightPlus;
  TH2D* hNUAweightMinus;
  // 15o
  TH3F* hCorrectNUAPos; // Protty
  TH3F* hCorrectNUANeg; // Protty
  ////////////////////////
  // VZERO
  ////////////////////////
  TList* fListVZEROCalib; // read list for V0 Calib
  // 10h
  TH2D* hMultV0Read;
  TProfile3D* pV0XMeanRead[3];
  TProfile3D* pV0YMeanRead[3];
  // 15o
  TH1D* hMultV0; // Dobrin
  AliOADBContainer* contMult;
  AliOADBContainer* contQxncm;
  AliOADBContainer* contQyncm;
  AliOADBContainer* contQxnam;
  AliOADBContainer* contQynam;
  TH1D* hQx2mV0[2];
  TH1D* hQy2mV0[2];
  TSpline3* splQ2c[90];

  // 18q/r
  TH2F* fHCorrectV0ChWeghts;
  TH2D* hQnPercentile;
  ////////////////////////
  // ZDC
  ////////////////////////
  TList* fListZDCCalib;
  // 10h
  TTree* tree;
  float vtxQuant1[3];
  float vtxQuant2[3];
  TProfile* fProfileForZNCGE;
  TProfile* fProfileForZNAGE;
  THnSparseF* fHn4DForZNCQxRC;
  THnSparseF* fHn4DForZNCQyRC;
  THnSparseF* fHn4DForZNCMtRC;
  THnSparseF* fHn4DForZNAQxRC;
  THnSparseF* fHn4DForZNAQyRC;
  THnSparseF* fHn4DForZNAMtRC;
  THnSparseI* fHn4DForZNCCountsRC;
  THnSparseI* fHn4DForZNACountsRC;
  TProfile2D* fProfile2DForCosC;
  TProfile2D* fProfile2DForSinC;
  TProfile2D* fProfile2DForCosA;
  TProfile2D* fProfile2DForSinA;
  // 15o 
  // 18q 18r
  TH1D* fHZDCCparameters;
  TH1D* fHZDCAparameters;

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
  // Psi QA
  // V0C [0]GE [1]RC
  TProfile* fProfileV0CQxCent[2];
  TProfile* fProfileV0CQyCent[2];
  TProfile* fProfileV0CQxVtx[2];
  TProfile* fProfileV0CQyVtx[2];
  TH2D* fHist2CalibPsi2V0CCent[2];
  TH3D* fHist3CalibQxQyCentV0C[2];
  TProfile* fProfileV0AQxCent[2];
  TProfile* fProfileV0AQyCent[2];
  TProfile* fProfileV0AQxVtx[2];
  TProfile* fProfileV0AQyVtx[2];
  TH2D* fHist2CalibPsi2V0ACent[2];
  TH3D* fHist3CalibQxQyCentV0A[2];
  // ZNC
  TProfile* fProfileZNCTowerMeanEnegry[2]; // ![0]Raw [1]GE
  TProfile* fProfileZNCQxCent[2];          // ![0]GE  [1]RC
  TProfile* fProfileZNCQyCent[2];          // ![0]GE  [1]RC
  TH2D* fHist2CalibPsi1ZNCCent[3];         // ![0]GE  [1]RC [2]SF
  // ZNA
  TProfile* fProfileZNATowerMeanEnegry[2];
  TProfile* fProfileZNAQxCent[2];
  TProfile* fProfileZNAQyCent[2];
  TH2D* fHist2CalibPsi1ZNACent[3];
  // ZNC-ZNA
  TProfile* fProfileZDCQxAQxCCent[2];
  TProfile* fProfileZDCQxAQyCCent[2];
  TProfile* fProfileZDCQyAQxCCent[2];
  TProfile* fProfileZDCQyAQyCCent[2];

  //Pion QA
  TH1D* fHistPionPt;
  TH1D* fHistPionEta;
  TH1D* fHistPionPhi;
  TH1D* fHistPionDcaXY;
  TH1D* fHistPionDcaZ;
  TH2D* fHist2PionSigTPC;
  TH2D* fHist2PionSigTOF;
  TH2D* fHist2PionSigRMS;
  TH3D* fHist3PionSigTPCTOFPt;
  TH1D* fHistAntiPionPt;
  TH1D* fHistAntiPionEta;
  TH1D* fHistAntiPionPhi;
  TH1D* fHistAntiPionDcaXY;
  TH1D* fHistAntiPionDcaZ;
  TH2D* fHist2AntiPionSigTPC;
  TH2D* fHist2AntiPionSigTOF;
  TH2D* fHist2AntiPionSigRMS;
  TH3D* fHist3AntiPionSigTPCTOFPt;

  //Kaon QA
  TH1D* fHistKaonPt;
  TH1D* fHistKaonEta;
  TH1D* fHistKaonPhi;
  TH1D* fHistKaonDcaXY;
  TH1D* fHistKaonDcaZ;
  TH2D* fHist2KaonSigTPC;
  TH2D* fHist2KaonSigTOF;
  TH2D* fHist2KaonSigRMS;
  TH3D* fHist3KaonSigTPCTOFPt;
  TH1D* fHistAntiKaonPt;
  TH1D* fHistAntiKaonEta;
  TH1D* fHistAntiKaonPhi;
  TH1D* fHistAntiKaonDcaXY;
  TH1D* fHistAntiKaonDcaZ;
  TH2D* fHist2AntiKaonSigTPC;
  TH2D* fHist2AntiKaonSigTOF;
  TH2D* fHist2AntiKaonSigRMS;
  TH3D* fHist3AntiKaonSigTPCTOFPt;

  //Proton QA
  TH1D* fHistProtonPt;
  TH1D* fHistProtonEta;
  TH1D* fHistProtonPhi;
  TH1D* fHistProtonDcaXY;
  TH1D* fHistProtonDcaZ;
  TH2D* fHist2ProtonSigTPC;
  TH2D* fHist2ProtonSigTOF;
  TH2D* fHist2ProtonSigRMS;
  TH3D* fHist3ProtonSigTPCTOFPt;
  TH1D* fHistAntiProtonPt;
  TH1D* fHistAntiProtonEta;
  TH1D* fHistAntiProtonPhi;
  TH1D* fHistAntiProtonDcaXY;
  TH1D* fHistAntiProtonDcaZ;
  TH2D* fHist2AntiProtonSigTPC;
  TH2D* fHist2AntiProtonSigTOF;
  TH2D* fHist2AntiProtonSigRMS;
  TH3D* fHist3AntiProtonSigTPCTOFPt;

  /////////////
  // Results //
  /////////////
  TList* fResultsList;
  // Plane
  TH2D* fHist2Psi2TPCPosCent;
  TH2D* fHist2Psi2TPCNegCent;
  TH2D* fHist2Psi2V0CCent;
  TH2D* fHist2Psi2V0ACent;
  TH2D* fHist2Psi1ZNCCent;
  TH2D* fHist2Psi1ZNACent;
  TH2D* fHist2Psi2TPCCent;

  // Res
  TProfile* fProfileTPCPsi2Correlation;       // TPCPos-TPCNeg
  TProfile* fProfileV0MPsi2Correlation;       // VOC-VOA
  TProfile* fProfileZDCPsi1Correlation;       // ZNC-ZNA 1st
  TProfile* fProfileZDCPsi2Correlation;       // ZNC-ZNA 2nd
  TProfile* fProfileV0CTPCPosPsi2Correlation; // V0C-TPCPos
  TProfile* fProfileV0ATPCPosPsi2Correlation; // V0A-TPCNeg
  TProfile* fProfileV0CTPCNegPsi2Correlation; // V0C-TPCNeg
  TProfile* fProfileV0ATPCNegPsi2Correlation; // V0A-TPCNeg
  TProfile* fProfileV0CTPCPsi2Correlation;    // V0C-TPC
  TProfile* fProfileV0ATPCPsi2Correlation;    // V0A-TPC

  // Flow
  // [0]TPC [1]V0C [2]V0A [3]ZNC [4]ZNA
  // [0]pos [1]neg(anti)
  TProfile2D* fProfile2RawFlowPtCentHadron[5][2];
  TProfile2D* fProfile2RawFlowPtCentProton[5][2];
  TProfile2D* fProfile2RawFlowPtCentPion[5][2];
  TProfile2D* fProfile2RawFlowPtCentKaon[5][2];
  TProfile2D* fProfile2RawFlowHadronQ2;
  TProfile2D* fProfile2RawFlowProtonQ2;
  TProfile2D* fProfile2RawFlowPionQ2;
  TProfile2D* fProfile2RawFlowKaonQ2;

  // δ
  TProfile* fProfileDeltaPionKaon[2]; //![0]:π-K same_charge  [1]:π-K opposit_charge
  TProfile* fProfileDeltaPionProton[2];
  TProfile* fProfileDeltaKaonProton[2];
  TProfile* fProfileDeltaPionPion[2];
  TProfile* fProfileDeltaKaonKaon[2];
  TProfile* fProfileDeltaProtonProton[2];
  TProfile* fProfileDeltaHadronHadron[2];
  // γ [5]:plane type [2]:pair type
  TProfile* fProfileGammaPionKaon[5][2];
  TProfile* fProfileGammaPionProton[5][2];
  TProfile* fProfileGammaKaonProton[5][2];
  TProfile* fProfileGammaPionPion[5][2];
  TProfile* fProfileGammaKaonKaon[5][2];
  TProfile* fProfileGammaProtonProton[5][2];
  TProfile* fProfileGammaHadronHadron[5][2];

  // for SS and OS cross check
  TProfile* fProfileDeltaPionKaonSplit[4]; //![0]:pi+ k+  [1]:π- K- [2]pi+ K- [3]pi- K+
  TProfile* fProfileDeltaPionProtonSplit[4];
  TProfile* fProfileDeltaKaonProtonSplit[4];
  TProfile* fProfileDeltaPionPionSplit[4];
  TProfile* fProfileDeltaKaonKaonSplit[4];
  TProfile* fProfileDeltaProtonProtonSplit[4];
  TProfile* fProfileDeltaHadronHadronSplit[4];
  // γ [5]:plane type [2]:pair type
  TProfile* fProfileGammaPionKaonSplit[5][4];
  TProfile* fProfileGammaPionProtonSplit[5][4];
  TProfile* fProfileGammaKaonProtonSplit[5][4];
  TProfile* fProfileGammaPionPionSplit[5][4];
  TProfile* fProfileGammaKaonKaonSplit[5][4];
  TProfile* fProfileGammaProtonProtonSplit[5][4];
  TProfile* fProfileGammaHadronHadronSplit[5][4];

  // Diff δ(ΔpT)
  TProfile2D* fProfile2DiffDeltaPionKaonDPt[2]; 
  TProfile2D* fProfile2DiffDeltaPionProtonDPt[2]; 
  TProfile2D* fProfile2DiffDeltaKaonProtonDPt[2]; 
  TProfile2D* fProfile2DiffDeltaPionPionDPt[2];
  TProfile2D* fProfile2DiffDeltaKaonKaonDPt[2]; 
  TProfile2D* fProfile2DiffDeltaProtonProtonDPt[2];
  TProfile2D* fProfile2DiffDeltaHadronHadronDPt[2];
  // Diff δ(SpT)
  TProfile2D* fProfile2DiffDeltaPionKaonSPt[2]; 
  TProfile2D* fProfile2DiffDeltaPionProtonSPt[2]; 
  TProfile2D* fProfile2DiffDeltaKaonProtonSPt[2]; 
  TProfile2D* fProfile2DiffDeltaPionPionSPt[2];
  TProfile2D* fProfile2DiffDeltaKaonKaonSPt[2]; 
  TProfile2D* fProfile2DiffDeltaProtonProtonSPt[2];
  TProfile2D* fProfile2DiffDeltaHadronHadronSPt[2];
  // Diff δ(Δη)
  TProfile2D* fProfile2DiffDeltaPionKaonDEta[2]; 
  TProfile2D* fProfile2DiffDeltaPionProtonDEta[2]; 
  TProfile2D* fProfile2DiffDeltaKaonProtonDEta[2]; 
  TProfile2D* fProfile2DiffDeltaPionPionDEta[2];
  TProfile2D* fProfile2DiffDeltaKaonKaonDEta[2]; 
  TProfile2D* fProfile2DiffDeltaProtonProtonDEta[2];
  TProfile2D* fProfile2DiffDeltaHadronHadronDEta[2];
  // Diff δ(q2);
  TProfile2D* fProfile2DiffDeltaPionKaonQ2[2]; 
  TProfile2D* fProfile2DiffDeltaPionProtonQ2[2]; 
  TProfile2D* fProfile2DiffDeltaKaonProtonQ2[2]; 
  TProfile2D* fProfile2DiffDeltaPionPionQ2[2];
  TProfile2D* fProfile2DiffDeltaKaonKaonQ2[2]; 
  TProfile2D* fProfile2DiffDeltaProtonProtonQ2[2];
  TProfile2D* fProfile2DiffDeltaHadronHadronQ2[2];

  // Diff γ(ΔpT)(only V0C Plane)
  TProfile2D* fProfile2DiffGammaPionKaonDPt[2]; 
  TProfile2D* fProfile2DiffGammaPionProtonDPt[2]; 
  TProfile2D* fProfile2DiffGammaKaonProtonDPt[2]; 
  TProfile2D* fProfile2DiffGammaPionPionDPt[2];
  TProfile2D* fProfile2DiffGammaKaonKaonDPt[2]; 
  TProfile2D* fProfile2DiffGammaProtonProtonDPt[2];
  TProfile2D* fProfile2DiffGammaHadronHadronDPt[2];
  // Diff γ(SpT)(only V0C Plane)
  TProfile2D* fProfile2DiffGammaPionKaonSPt[2]; 
  TProfile2D* fProfile2DiffGammaPionProtonSPt[2]; 
  TProfile2D* fProfile2DiffGammaKaonProtonSPt[2]; 
  TProfile2D* fProfile2DiffGammaPionPionSPt[2];
  TProfile2D* fProfile2DiffGammaKaonKaonSPt[2]; 
  TProfile2D* fProfile2DiffGammaProtonProtonSPt[2];
  TProfile2D* fProfile2DiffGammaHadronHadronSPt[2];
  // Diff γ(Δη)
  TProfile2D* fProfile2DiffGammaPionKaonDEta[2]; 
  TProfile2D* fProfile2DiffGammaPionProtonDEta[2]; 
  TProfile2D* fProfile2DiffGammaKaonProtonDEta[2]; 
  TProfile2D* fProfile2DiffGammaPionPionDEta[2];
  TProfile2D* fProfile2DiffGammaKaonKaonDEta[2]; 
  TProfile2D* fProfile2DiffGammaProtonProtonDEta[2];
  TProfile2D* fProfile2DiffGammaHadronHadronDEta[2];
  // Diff γ(q2)
  TProfile2D* fProfile2DiffGammaPionKaonQ2[2]; 
  TProfile2D* fProfile2DiffGammaPionProtonQ2[2]; 
  TProfile2D* fProfile2DiffGammaKaonProtonQ2[2]; 
  TProfile2D* fProfile2DiffGammaPionPionQ2[2];
  TProfile2D* fProfile2DiffGammaKaonKaonQ2[2]; 
  TProfile2D* fProfile2DiffGammaProtonProtonQ2[2];
  TProfile2D* fProfile2DiffGammaHadronHadronQ2[2];


  // C(Δη,Δφ) [cent][pair type]
  TH2D* fHist2DEtaDPhiPionKaon[8][2];
  TH2D* fHist2DEtaDPhiPionProton[8][2];
  TH2D* fHist2DEtaDPhiKaonProton[8][2];
  TH2D* fHist2DEtaDPhiPionPion[8][2];
  TH2D* fHist2DEtaDPhiKaonKaon[8][2];
  TH2D* fHist2DEtaDPhiProtonProton[8][2];
  TH2D* fHist2DEtaDPhiHadronHadron[8][2];

  // C(Δη,sφ) [cent][pair type] only V0C Plane
  TH2D* fHist2DEtaSPhiPionKaon[8][2];
  TH2D* fHist2DEtaSPhiPionProton[8][2];
  TH2D* fHist2DEtaSPhiKaonProton[8][2];
  TH2D* fHist2DEtaSPhiPionPion[8][2];
  TH2D* fHist2DEtaSPhiKaonKaon[8][2];
  TH2D* fHist2DEtaSPhiProtonProton[8][2];
  TH2D* fHist2DEtaSPhiHadronHadron[8][2];

  AliAnalysisTaskPIDCME(const AliAnalysisTaskPIDCME&);
  AliAnalysisTaskPIDCME& operator=(const AliAnalysisTaskPIDCME&);

  ClassDef(AliAnalysisTaskPIDCME, 1);
};

#endif