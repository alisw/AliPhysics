#ifndef ALIANALYSISTASKITSSATRACKS
#define ALIANALYSISTASKITSSATRACKS

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskITSsaTracks
// AliAnalysisTaskSE to extract QA and performance histos for ITS standalone tracks
// 
//
// Authors: L. Milano, milano@to.infn.it
//          F. Prino, prino@to.infn.it
//          
//*************************************************************************

class TList;
class TH1F;
class TH2F;
class TTree;
class TString;
class AliESDEvent;
class AliESDfriend;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskITSsaTracks : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskITSsaTracks();
  virtual ~AliAnalysisTaskITSsaTracks();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

  void SetMinITSPoints(Int_t minp=4){
    fMinITSpts=minp;
  }
  void SetMinTPCPoints(Int_t minp=50){
    fMinTPCpts=minp;
  }
  void SetMinSPDPoints(Int_t minp=1){
    fMinSPDpts=minp;
  }
  void SetMinPointsForITSPid(Int_t minp=3){
    fMinITSptsForPid=minp;
  }
  void SetMinPointsForMatchToTPC(Int_t minp=6){
    fMinITSptsForMatch=minp;
  }

  void SetPtBins(Int_t n, Double_t* lim);

  void RequirePointInLayer(Int_t iLay){
    if(iLay<6) fRequirePoint[iLay]=kTRUE;
  }
  void RequireInnerSPD(){
    fRequirePoint[0]=kTRUE;
  }
  void RequireBothSPD(){
    fRequirePoint[0]=kTRUE;
    fRequirePoint[1]=kTRUE;
  }

  void SetReadMC(Bool_t optMC=kTRUE){
    fReadMC=optMC;
  }
  void SetUseMCtruthForPID(Bool_t opt=kTRUE){
    fUseMCId=opt;
  }


 private:
  enum {kPion=0,kKaon,kProton,kNspecies};
  enum {kMaxPtBins=40};

  AliAnalysisTaskITSsaTracks(const AliAnalysisTaskITSsaTracks &source);
  AliAnalysisTaskITSsaTracks& operator=(const AliAnalysisTaskITSsaTracks &source);
  
  TList*  fOutput;          //! list of output histos
  TH1F*   fHistNEvents;     //! histo with N of events  

  TH1F*   fHistPtTPCITS[kNspecies];    //! pt distribution of TPC+ITS tracks
  TH1F*   fHistPtITSsa[kNspecies];     //! pt distribution of ITSsa tracks
  TH1F*   fHistPtITSpureSA[kNspecies]; //! pt distribution of ITS pure SA tracks

  TH1F*   fHistPtTPCITSAll;     //! pt distribution of all TPC+ITS tracks 
  TH1F*   fHistPtTPCITSGood;    //! pt distribution of good TPC+ITS tracks 
  TH1F*   fHistPtTPCITSFake;    //! pt distribution of fake TPC+ITS tracks
  TH1F*   fHistPtITSsaAll;      //! pt distribution of all ITSsa tracks
  TH1F*   fHistPtITSsaGood;     //! pt distribution of good ITSsa tracks
  TH1F*   fHistPtITSsaFake;     //! pt distribution of fake ITSsa tracks
  TH1F*   fHistPtITSpureSAAll;  //! pt distribution of all ITS pure SA tracks
  TH1F*   fHistPtITSpureSAGood; //! pt distribution of good ITS pure SA tracks
  TH1F*   fHistPtITSpureSAFake; //! pt distribution of fake ITS pure SA tracks

  TH2F*   fHistEtaPhiTPCITS[kNspecies];    //! etaphi distr. of TPC+ITS tracks
  TH2F*   fHistEtaPhiITSsa[kNspecies];     //! etaphi distr. of ITSsa tracks
  TH2F*   fHistEtaPhiITSpureSA[kNspecies]; //! etaphi distr. of ITSpureSA tracks

  TH2F*   fHistdedxvsPtITSpureSA3cls; //! dedx for ITSpureSA tracks vs. pt
  TH2F*   fHistdedxvsPITSpureSA3cls; //! dedx for ITSpureSA tracks vs. p
  TH2F*   fHistdedxvsPtITSpureSA4cls; //! dedx for ITSpureSA tracks vs. pt
  TH2F*   fHistdedxvsPITSpureSA4cls; //! dedx for ITSpureSA tracks vs. p

  TH2F*   fHistNcluTPCITS[kNspecies];    //! n. of clusters for TPC+ITS tracks vs. pt
  TH2F*   fHistNcluITSsa[kNspecies];     //! n. of clusters for ITSsa tracks vs. pt
  TH2F*   fHistNcluITSpureSA[kNspecies]; //! n. of clusters for ITSpureSA tracks vs. pt
  TH2F*   fHistd0rphiITSpureSA[kNspecies]; //! d0z for ITSpureSA tracks vs. pt
  TH2F*   fHistd0zITSpureSA[kNspecies]; //! d0z for ITSpureSA tracks vs. pt
  TH2F*   fHistCluInLayTPCITS[kNspecies];    //! TPC+ITS tracks with cluster in layer 
  TH2F*   fHistCluInLayITSsa[kNspecies];     //! ITSsa tracks with cluster in layer 
  TH2F*   fHistCluInLayITSpureSA[kNspecies]; //! for ITSpureSA tracks with cluster in layer

  TH2F*   fHistOuterLayITSpureSA[kNspecies]; //! outer layer with cluster vs. pt

  TH2F*   fHistPtResid[kNspecies];         //! pt residuals (TPC) in pt bins
  TH2F*   fHistPtRelResid[kNspecies];      //! pt relative residuals (TPC) in pt bins
  TH2F*   fHistInvPtResid[kNspecies];      //! 1/pt residuals (TPC) in pt bins
  TH2F*   fHistInvPtRelResid[kNspecies];   //! 1/pt relative resid. (TPC) in pt bins
  TH2F*   fHistMCPtResid[kNspecies];       //! pt residuals (MC) vs. pt     
  TH2F*   fHistMCPtRelResid[kNspecies];    //! pt relative residuals (MC) vs. pt     
  TH2F*   fHistMCInvPtResid[kNspecies];    //! 1/pt residuals (MC) vs. pt
  TH2F*   fHistMCInvPtRelResid[kNspecies]; //! 1/pt relative residulas (MC) vs. pt

  Int_t fNPtBins;                  // number of Pt bins
  Float_t fPtLimits[kMaxPtBins+1]; // Pt bin limits
  Int_t   fMinITSpts;          // Minimum number of ITS points per track
  Int_t   fMinSPDpts;          // Minimum number of SPD points per track
  Int_t   fMinITSptsForPid;    // Minimum number of SDD+SSD points per track
  Int_t   fMinITSptsForMatch;  // Minimum number of ITS points to macth to TPC
  Int_t   fMinTPCpts;          // Minimum number of TPC points per track
  Bool_t  fRequirePoint[6]; // require point in given layer
  Bool_t  fReadMC;          // flag read/not-read MC truth info
  Bool_t  fUseMCId;         // flag use/not-use MC identity for PID

  ClassDef(AliAnalysisTaskITSsaTracks,1);  
};


#endif

