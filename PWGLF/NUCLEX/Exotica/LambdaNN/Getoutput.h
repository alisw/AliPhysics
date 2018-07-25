#ifndef Getoutput_H
#define Getoutput_H

#include <TList.h>
#include <TFile.h>
#include <TSystem.h>
#include <TNtupleD.h>
#include <TObjString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TDatabasePDG.h>
#include <TMath.h>

//using namespace std;

class Getoutput {
 public:
  Getoutput();
  bool LoadParams(const char *paramfile="params.txt");
  bool LoadFile(const char *inputfile="AnalysisResults.root",const char *listName="clistNtupHyper");
  void LoopOverV0(Int_t Hcharge=-1);
  void StoreOutputData(const char *filename="results.root");
  bool LoadOutputData(const char *filename="results3H1.root");
   void DrawResults();
  void ClearInputData();

  void SetMC(){fIsMC=kTRUE;}
  void SetAODCuts(){fUseAODCut=kTRUE;}
  void SetTOFpid(){fIncludePidTOF=kTRUE;}
  void Set3HPcut(Double_t pMin){f3HPcut=pMin;}
  void BookOutputData();
  Double_t GetInvMass (TVector3 vPos, TVector3 vNeg, Double_t mPos, Double_t mNeg);
  bool EventSelectionAOD(Double_t *arr);

  Bool_t fIsMC;
  Bool_t fIncludePidTOF;
  Bool_t fRejectBkg; // useful in case like-sign V0 are produced
  Bool_t fUseAODCut;
  Int_t f3Hsign;
  Double_t f3HPcut;
  Double_t array[28];
  Float_t param[20];

  TList *fInputList;
  TList *fHistList;
  TList *fOutputList;

  TH2F *hDiffRelP[6];
  TH1F *hDecayLength[5];
  TH1F *hProdVtx[5];
  TH1F *hMass[5];
  TH1F *hMassTrd[5];
  TH1F *hMassBkg;
  TH1F *hMassSignal;
  TH1F *hMassContrib[2][7]; // quark u&d, s,c,b for pion and triton
  TH1F *hMassK0;
  TH1F *hMassLambda;
  TH1F *hMassGamma;
  TH1F *hHighM;
  TH1F *hTriP[5];
  TH1F *hPiP[5];
  TH1F *hV0P[5];
  TH2F *triTOFmass;
  TH1F *hDcaD;
  TH2F *hArmPlot;
  TH2F *hArmPlotSel[5];
  TH2F *hTPCsignalPi;
  TH2F *hTPCsignalTri;
  TH2F *hTPCsignalTriAll;
  TH2F *hTPCsignalPiClean;
  TH2F *hTPCsignalTriClean;
  TH2F* hTPCsignalTri91Lim;
  TH2F* hTPCsignalTriTrd;
  TH2F *hMumCheck[2];
  TH1I *hMonitorPlot;
 

  TNtupleD *ntTot;

  static const Int_t fgArrSize;
  static Double_t fgLH3Mass;
  static Double_t fgPionMass;
  static Double_t fgEleMass;
  static Double_t fgProtMass;

};


enum { kPposx, kPposy, kPposz, kPnegx, kPnegy, kPnegz,	//0-5
 kNSPi, kNSTri, kTriTOFmass, kPiTPCsignal, kTriTPCsignal,	// 6-10
 kV0mom, kPtArm, kAlphaArm,	// 11-13
 kDcaTriXY, kDcaTriZ, kV0dcaD, kDecayPath, kDecayPathXY,	// 14-18
 kV0Dca, kCosP, kV0VtxErrSum, kSign,	// 19-22
 kDcaPi, kIsTrdEle, kSigPiFromPiTof, kSigPrTof, kSigPiTof, kNclusITS,
 kPiPdgCode, kTriPdgCode, kMumPiPdgCode, kMumTriPdgCode
};				//23-27

enum {
 kParMinP, kParMinPv0, kParMaxP3H, kParPiLim, kPar3hLim, kParNclusITS, kParNsigmaPID, kParNsigmaTOFmass, kParDcaTriZ, kParCosP, kParV0Dca,
 kK0MassLow, kK0MassHigh, kLambdaMassLow, kLambdaMassHigh,kGammaMassHigh,kTOFpid,kIsMc,kPIDResponseYear,k3HPlim
};
#endif

   
