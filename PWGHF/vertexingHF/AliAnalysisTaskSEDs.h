#ifndef ALIANALYSISTASKDS_H
#define ALIANALYSISTASKDS_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
///                                                                      //
/// \class Analysis task to produce Ds candidates mass spectra           //
/// \author Origin: F.Prino, Torino, prino@to.infn.it                    //
///                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliLog.h"

class AliNormalizationCounter;


class AliAnalysisTaskSEDs : public AliAnalysisTaskSE
{
 public:
 
  AliAnalysisTaskSEDs();
  AliAnalysisTaskSEDs(const char *name, AliRDHFCutsDstoKKpi* analysiscuts, Int_t fillNtuple=0);
  virtual ~AliAnalysisTaskSEDs();
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}	
  void SetWriteOnlySignalInNtuple(Bool_t opt=kTRUE){
    if(fReadMC) fWriteOnlySignal=opt;
    else AliError("fReadMC has to be kTRUE");
  }
  void SetFillNtuple(Int_t fill=0){fFillNtuple=fill;}
  void SetFillNSparse(Bool_t fill=kTRUE){fFillSparse=fill;}
  void SetMassRange(Double_t rang=0.4){fMassRange=rang;}
  void SetDoCutVarHistos(Bool_t opt=kTRUE) {fDoCutVarHistos=opt;}
  void SetUseSelectionBit(Bool_t opt=kFALSE){ fUseSelectionBit=opt;}
  void SetAODMismatchProtection(Bool_t opt=kTRUE) {fAODProtection=opt;}
  Bool_t CheckDaugAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
  void FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader);
  
  void SetInvMassBinSize(Double_t binsiz=0.002){fMassBinSize=binsiz;}
  void SetPtBins(Int_t n, Float_t* lim);
  void SetAnalysisCuts(AliRDHFCutsDstoKKpi* cuts){fAnalysisCuts=cuts;}
  void SetSystem(Int_t system){fSystem = system;}
  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*4;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*4+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*4+2;}
  Int_t GetReflSignalHistoIndex(Int_t iPtBin) const { return iPtBin*4+3;}

  enum {kMaxPtBins=20,knVarForSparse=13,knVarForSparseAcc=2,knVarForSparseIP=6};

  AliAnalysisTaskSEDs(const AliAnalysisTaskSEDs &source);
  AliAnalysisTaskSEDs& operator=(const AliAnalysisTaskSEDs& source); 

  TList*  fOutput;                    //!<! list send on output slot 0
  TH1F*   fHistNEvents;               //!<! hist. for No. of events  
  TH1F*   fChanHist[4];               //!<! hist. with KKpi and piKK candidates (sig,bkg,tot)
  TH1F*   fMassHist[4*kMaxPtBins];    //!<! hist. of mass spectra (sig,bkg,tot)
  TH1F*   fMassHistPhi[4*kMaxPtBins];    //!<! hist. of mass spectra via phi (sig,bkg,tot)
  TH1F*   fMassHistK0st[4*kMaxPtBins];    //!<! hist. of mass spectra via K0* (sig,bkg,tot)
  TH1F*   fMassHistKK[kMaxPtBins];    //!<! hist. of mass spectra of KK
  TH1F*   fMassHistKpi[kMaxPtBins];    //!<! hist. of mass spectra of Kpi
  TH1F*   fCosPHist[4*kMaxPtBins];    //!<! hist. of cos pointing angle (sig,bkg,tot)
  TH1F*   fDLenHist[4*kMaxPtBins];    //!<! hist. of decay length (sig,bkg,tot)
  TH1F*   fSumd02Hist[4*kMaxPtBins];  //!<! hist. for sum d02 (Prod Cuts)
  TH1F*   fSigVertHist[4*kMaxPtBins]; //!<! hist. for sigVert (Prod Cuts)
  TH1F*   fPtMaxHist[4*kMaxPtBins];   //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fPtCandHist[4*kMaxPtBins];  //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fDCAHist[4*kMaxPtBins];     //!<! hist. for DCA (Prod Cuts)
  TH1F*   fPtProng0Hist[4*kMaxPtBins]; //!<! hist. for Pt Max (Prod Cuts)
  TH1F*   fPtProng1Hist[4*kMaxPtBins]; //!<! hist. for DCA (Prod Cuts)
  TH1F*   fPtProng2Hist[4*kMaxPtBins]; //!<! hist. for DCA (Prod Cuts)
  TH2F*   fDalitz[4*kMaxPtBins];      //!<! dalitz plot (sig,bkg,tot)
  TH2F*   fDalitzPhi[4*kMaxPtBins];   //!<! dalitz plot via phi (sig,bkg,tot)
  TH2F*   fDalitzK0st[4*kMaxPtBins];   //!<! dalitz plot via K0* (sig,bkg,tot)
  TH2F *fPtVsMass;    //!<! hist. of pt vs. mass (prod. cuts)
  TH2F *fPtVsMassPhi;    //!<! hist. of pt vs. mass (phi selection)
  TH2F *fPtVsMassK0st;   //!<! hist. of pt vs. mass (K0* selection)
  TH2F *fYVsPt;       //!<! hist. of Y vs. Pt (prod. cuts)
  TH2F *fYVsPtSig;    //!<! hist. of Y vs. Pt (MC, only sig, prod. cuts)
  TH1F *fHistCentrality[3];//!<!hist. for cent distr (all,sel ev, )
  TH2F *fHistCentralityMult[3];//!<!hist. for cent distr vs mult (all,sel ev, )
  TNtuple *fNtupleDs; //!<! output ntuple
  Int_t fFillNtuple;                 /// 0 not to fill ntuple
                                     /// 1 for filling ntuple for events through Phi
                                     /// 2 for filling ntuple for events through K0Star
                                     /// 3 for filling all
  Int_t   fSystem;                    /// 0 = pp, 1 = pPb,PbPb
  Bool_t  fReadMC;                    ///  flag for access to MC
  Bool_t  fWriteOnlySignal;           ///  flag to control ntuple writing in MC
  Bool_t  fDoCutVarHistos;            ///  flag to create and fill histos with distributions of cut variables
  Bool_t  fUseSelectionBit;           /// flag for usage of HasSelectionBit
  Bool_t  fFillSparse;                /// flag for usage of THnSparse
  Bool_t  fAODProtection;             /// flag to activate protection against AOD-dAOD mismatch
  UChar_t fNPtBins;                   /// number of Pt bins
  TList *fListCuts; //list of cuts
  Float_t fPtLimits[kMaxPtBins+1];    ///  limits for pt bins
  Double_t fMassRange;                /// range for mass histogram
  Double_t fMassBinSize;              /// bin size for inv. mass histo

  AliNormalizationCounter *fCounter;//!<!Counter for normalization
  AliRDHFCutsDstoKKpi *fAnalysisCuts; /// Cuts for Analysis
  
  THnSparseF *fnSparse;       ///!<!THnSparse for candidates on data
  THnSparseF *fnSparseIP;       ///!<!THnSparse for topomatic variable
  THnSparseF *fnSparseMC[4];  ///!<!THnSparse for MC
  ///[0]: Acc step prompt Ds
  ///[1]: Acc step FD Ds
  ///[2]: Selected prompt Ds
  ///[3]: Selected FD Ds
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDs,17);    ///  AliAnalysisTaskSE for Ds mass spectra
  /// \endcond
};

#endif


