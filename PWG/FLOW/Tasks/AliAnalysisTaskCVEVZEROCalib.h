#ifndef AliAnalysisTaskCVEVZEROCalib_cxx
#define AliAnalysisTaskCVEVZEROCalib_cxx
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


class AliAnalysisTaskCVEVZEROCalib : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskCVEVZEROCalib();
  AliAnalysisTaskCVEVZEROCalib(const char* name);
  virtual ~AliAnalysisTaskCVEVZEROCalib();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  void IfDebug(bool bDebug) { this->fDebug = bDebug; }
  // read in
  void SetListForVZEROCalib(TList* flist) { this->fListVZEROCalib = (TList*)flist->Clone(); }

  // Global
  void SetTrigger(TString trigger) { this->fTrigger = trigger; }
  void SetPeriod(TString period) { this->fPeriod = period; }
  // Event
  void SetVzCut(double vzCut) { this->fVzCut = vzCut; }
  void SetCentCut(float centDiffCut) { this->fCentDiffCut = centDiffCut; }

  void SetVzBinNumber(int num) {this->fVzBinNum = num; }
  void SetCentBinNumber(int num) {this->fCentBinNum = num; }

 private:

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
  bool RecenterVZEROQVector();

  //////////////////////
  // Cuts and options //
  //////////////////////
  // Global
  TString fTrigger; //
  TString fPeriod;  // period
  // Event
  double fVzCut;      // vz cut
  float fCentDiffCut; // centrality restriction for V0M and TRK

  ///////////////////The following files are from the data//////////////////////////////////
  /////////////
  // Handles //
  /////////////
  AliAODEvent* fAOD;            // aod Event
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
  int fVzBinNum;
  int fCentBinNum;

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
  // 18q
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
  TProfile* fProfileChanalMult[2];
  TProfile* fProfileQ2xV0CCent[2];
  TProfile* fProfileQ2yV0CCent[2];
  TProfile* fProfileQ2xV0ACent[2];
  TProfile* fProfileQ2yV0ACent[2];
  TProfile* fProfileQ2xV0CVz[2];
  TProfile* fProfileQ2yV0CVz[2];
  TProfile* fProfileQ2xV0AVz[2];
  TProfile* fProfileQ2yV0AVz[2];

  /////////////
  // Results //
  /////////////
  TList* fResultsList;
  TProfile2D* p2D_V0C_meanQ2x_cent_vz[150];
  TProfile2D* p2D_V0C_meanQ2y_cent_vz[150];
  TProfile2D* p2D_V0A_meanQ2x_cent_vz[150];
  TProfile2D* p2D_V0A_meanQ2y_cent_vz[150];

  AliAnalysisTaskCVEVZEROCalib(const AliAnalysisTaskCVEVZEROCalib&);
  AliAnalysisTaskCVEVZEROCalib& operator=(const AliAnalysisTaskCVEVZEROCalib&);

  ClassDef(AliAnalysisTaskCVEVZEROCalib, 1);
};

#endif

