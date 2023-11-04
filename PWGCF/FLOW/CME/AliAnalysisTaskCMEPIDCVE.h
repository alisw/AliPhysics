#ifndef AliAnalysisTaskCMEPIDCVE_cxx
#define AliAnalysisTaskCMEPIDCVE_cxx
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


class AliAnalysisTaskCMEPIDCVE : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskCMEPIDCVE();
  AliAnalysisTaskCMEPIDCVE(const char* name);
  virtual ~AliAnalysisTaskCMEPIDCVE();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  // Switch
  void IfDebug(bool bDebug) { this->fDebug = bDebug; }
  void IfTightPileUp(bool bTightPileUp) { this->isTightPileUp = bTightPileUp; }

  // Global
  void SetTrigger(TString trigger) { this->fTrigger = trigger; }
  void SetPeriod(TString period) { this->fPeriod = period; }

  // Event
  void SetVzCut(double vzCut) { this->fVzCut = vzCut; }
  void SetCentCut(float centDiffCut) { this->fCentDiffCut = centDiffCut; }

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

 private:
  ////////////////////////
  // Procedural function
  ////////////////////////
  void ResetVectors();
  bool LoopTracks();
  bool PairTrkTrk();

  ////////////////////////
  // Functional function
  ////////////////////////
  // Read in
  // Pile-up
  bool RejectEvtTFFit();
  // Track
  bool AcceptAODTrack(AliAODTrack* track);
  int GetPIDofParticle(AliAODTrack* ftrack);
  // Get DCA
  bool GetDCA(double &dcaxy, double &dcaz, AliAODTrack* track);

  ////////////////////////
  // input variables
  ////////////////////////
  // Switch
  bool isTightPileUp;

  // Global
  TString fTrigger;
  TString fPeriod;

  // Event
  double fVzCut;      // vz cut
  float fCentDiffCut; // centrality restriction for V0M and TRK

  // Track
  int fFilterBit;   // AOD filter bit selection
  int fNclsCut;     // ncls cut for all tracks
  float fChi2Max;   // upper limmit for chi2
  float fChi2Min;   // lower limmit for chi2
  float fDcaCutXY;  // dcaxy cut for all tracks
  float fDcaCutZ;   // dcaz cut for all tracks
  float fPtMin;     // minimum pt for tracks
  float fPtMax;     // maximum pt for tracks
  float fEtaCut;    // eta cut for tracks
  float fDedxCut;   // dedx cut for tracks

  // PID
  float fNSigmaTPCCut;
  float fNSigmaTOFCut;

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
  float fCent;       // centrality
  int fCentBin;      // centrality bin: 0-10
  float fCentV0M;    // Centrality V0M
  float fCentTRK;    // Centrality TRK
  float fCentSPD0;   // Centrality SPD0
  float fCentSPD1;   // Centrality SPD1

  // Vector for particles from Tracks [pt,eta,phi,pid,label]
  std::vector<std::array<double,5>> vecParticle;

  ///////////////////The following files are read from external sources////////////////////
  ////////////////////////
  // Pile up Function
  ////////////////////////
  TF1* fSPDCutPU;
  TF1* fV0CutPU;
  TF1* fCenCutLowPU;
  TF1* fCenCutHighPU;
  TF1* fMultCutPU;

  ///////////////////The following files will be saved//////////////////////////////////
  //////////////
  // QA Plots //
  //////////////
  TList* listQA;
  // General QA
  // Event-wise
  TH1D* hEvtCount;
  TH1I* hRunNumBin;

  TH1D* hCent[2];
  TH1D* hVz[2];
  TH2D* h2Cent[8];
  TH2D* h2MultCent[2];
  TH2D* h2MultMult[6];
  // Track-wise
  TH1D* hPt;
  TH1D* hEta;
  TH1D* hPhi;
  TH1D* hNhits;
  TH2D* h2PDedx;
  TH1D* hDcaXY;
  TH1D* hDcaZ;

  //[0]:matter [1]:antimatter
  //pion
  TH1D* hPtPion[2];
  TH1D* hEtaPion[2];
  TH1D* hPhiPion[2];
  //kaon
  TH1D* hPtKaon[2];
  TH1D* hEtaKaon[2];
  TH1D* hPhiKaon[2];
  //proton
  TH1D* hPtProton[2];
  TH1D* hEtaProton[2];
  TH1D* hPhiProton[2];

  /////////////
  // Results //
  /////////////
  TList* listResults;
  //        pion kaon proton
  //pion    00   01   02
  //kaon    10   11   12
  //proton  20   21   22

  TProfile* pDeltaOS[3][3];
  TProfile* pDeltaSS[3][3];

  AliAnalysisTaskCMEPIDCVE(const AliAnalysisTaskCMEPIDCVE&);
  AliAnalysisTaskCMEPIDCVE& operator=(const AliAnalysisTaskCMEPIDCVE&);

  ClassDef(AliAnalysisTaskCMEPIDCVE, 1);
};

#endif

