#ifndef ALIANALYSISTASKCHECKHFMCPROD
#define ALIANALYSISTASKCHECKHFMCPROD

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class Class AliAnalysisTaskCheckHFMCProd
/// \brief AliAnalysisTask to check MC production at ESD+Kine level
///
///
/// \author Author: F. Prino, prino@to.infn.it
///
//*************************************************************************

class TList;
class TNtuple;
class TH1F;
class TH2F;
class TH3F;
class TTree;
class TString;
class AliESDEvent;
class AliESDfriend;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckHFMCProd : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskCheckHFMCProd();
  virtual ~AliAnalysisTaskCheckHFMCProd();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

  void SetSearchUpToQuark(Bool_t opt){fSearchUpToQuark=opt;};


  void SetReadMC(Bool_t opt) {fReadMC=opt;}
  void SetpPb() {fSystem=2;}
  void SetPbPb() {fSystem=1;}
  void Setpp() {fSystem=0;}

 private:

  AliAnalysisTaskCheckHFMCProd(const AliAnalysisTaskCheckHFMCProd &source);
  AliAnalysisTaskCheckHFMCProd& operator=(const AliAnalysisTaskCheckHFMCProd &source);
  
  TList* fOutput;          //!<! list of output histos
  TH1F* fHistoNEvents;     //!<! histo with N of events

  TH1F* fHistoPhysPrim;     //!<! histo of n. of physical primaries in |eta|<0.5
  TH1F* fHistoPhysPrimPiKPi09; //!<! histo of n. of primary pi, K, p
  TH2F* fHistoPhysPrimPiKPi09vsb; //!<! histo of n. of primary pi, K, p vs. b
  TH1F* fHistoTracks;       //!<! histo with number of ESD tracks
  TH1F* fHistoSelTracks;    //!<! histo with number of TPC+ITS refit ESD tracks
  TH1F* fHistoTracklets;    //!<! histo with number of SPD tracklets
  TH1F* fHistoTrackletsEta1;//!<! histo with number of SPD tracklets in |eta|<1
  TH1F* fHistoPtPhysPrim;   //!<! histo of pt distribution of phys primaries
  TH1F* fHistoEtaPhysPrim;  //!<! histo of eta distribution of physical primaries
  
  TH1F* fHistoSPD3DVtxX;     //!<! histo with distr. of x coord. of SPD3D vertex
  TH1F* fHistoSPD3DVtxY;     //!<! histo with distr. of y coord. of SPD3D vertex
  TH1F* fHistoSPD3DVtxZ;     //!<! histo with distr. of z coord. of SPD3D vertex

  TH1F* fHistoSPDZVtxX;     //!<! histo with distr. of x coord. of SPDZ vertex
  TH1F* fHistoSPDZVtxY;     //!<! histo with distr. of y coord. of SPDZ vertex
  TH1F* fHistoSPDZVtxZ;     //!<! histo with distr. of z coord. of SPDZ vertex

  TH1F* fHistoTRKVtxX;     //!<! histo with distr. of x coord. of TRK vertex
  TH1F* fHistoTRKVtxY;     //!<! histo with distr. of y coord. of TRK vertex
  TH1F* fHistoTRKVtxZ;     //!<! histo with distr. of z coord. of TRK vertex

  TH2F* fHistoNcharmed;   //!<! histo of D mesons vs. dN/dy
  TH2F* fHistoNbVsNc;     //!<! histo of n. b quarks vs. n c. quarks

  TH2F*  fHistBYPtAllDecay[5];   //!<! histo of y vs. pt from prompt B0, B+, B*, Bs, Lb
  TH2F*  fHistYPtAllDecay[5];   //!<! histo of y vs. pt from prompt D0, D+, D*, Ds, Lc, no selection on decay channel  
  TH2F*  fHistYPtPromptAllDecay[5];   //!<! histo of y vs. pt from prompt D0, D+, D*, Ds, Lc, no selection on decay channel  
  TH2F*  fHistYPtFeeddownAllDecay[5];   //!<! histo of y vs. pt from prompt D0, D+, D*, Ds, Lc, no selection on decay channel
  TH2F*  fHistYPtPrompt[5];   //!<! histo of y vs. pt from prompt D0, D+, D*, Ds, Lc
  TH2F*  fHistYPtFeeddown[5]; //!<! histo of y vs. pt from feeddown D0, D+, D*, Ds, Lc
  TH2F* fHistYPtD0byDecChannel[2]; //!<! histo of y vs. pt for D0->Kpi and D0->Kpipipi
  TH2F* fHistYPtDplusbyDecChannel[3]; //!<! histo of y vs. pt for D+->Kpipi and D+->K0*pi
  TH2F* fHistYPtDsbyDecChannel[2]; //!<! histo of y vs. pt for Ds->phipi and Ds->K0*K
  TH1F* fHistOriginPrompt;    //!<! histo of D production point (prompt)
  TH1F* fHistOriginFeeddown;  //!<! histo of D production point (feeddown)
  TH1F* fHistMotherID;        //!<! histo of mother ID
  TH1F* fHistDSpecies;          //!<! histo of D hadron species
  TH1F* fHistBSpecies;          //!<! histo of B hadron species
  TH1F* fHistLcDecayChan;      //!<! histo of Lc decay modes
  TH2F* fHistNcollHFtype;      //!<! histo of number of injected events vs. type
  TH2F* fHistNinjectedvsb;     //!<! histo of number of injected events vs. b
  TH3F* fHistEtaPhiPtGenEle;   //!<! histo of generated electrons
  TH3F* fHistEtaPhiPtGenPi;   //!<! histo of generated pions
  TH3F* fHistEtaPhiPtGenK;   //!<! histo of generated kaons
  TH3F* fHistEtaPhiPtGenPro;   //!<! histo of generated protons
  TH3F* fHistEtaPhiPtRecEle;   //!<! histo of generated electrons
  TH3F* fHistEtaPhiPtRecPi;   //!<! histo of generated pions
  TH3F* fHistEtaPhiPtRecK;   //!<! histo of generated kaons
  TH3F* fHistEtaPhiPtRecPro;   //!<! histo of generated protons
  TH2F* fHistPtRecVsPtGen;   //!<! correlation between rec and gen pt
  TH2F* fHistPhiRecVsPhiGen;   //!<! correlation between rec and gen pt
  TH2F* fHistEtaRecVsEtaGen;   //!<! correlation between rec and gen pt
  TH1F* fHistPtRecGood;   //!<! pt distribution of "good" tracks
  TH1F* fHistPtRecFake;   //!<! pt distribution of fake tracks
  Bool_t fSearchUpToQuark; /// c/b separation using quarks
  Int_t fSystem;         /// 0=pp, 1=PbPb, 2=pPb
  AliESDtrackCuts *fESDtrackCuts; /// track selection
  Bool_t fReadMC;
    
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCheckHFMCProd,9);
  /// \endcond
};


#endif
