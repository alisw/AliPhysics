#ifndef ALIANALYSISBTOJPSITOELE_H
#define ALIANALYSISBTOJPSITOELE_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                         Class AliAnalysisBtoJPSItoEle
//                  Unbinned log-likelihood fit analysis class
//
//                             Origin: C.Di Giglio
//       Contact: Carmelo.Digiglio@ba.infn.it , giuseppe.bruno@ba.infn.it
//-------------------------------------------------------------------------
class TNtuple ;
class AliBtoJPSItoEleCDFfitFCN ;
class AliBtoJPSItoEleCDFfitHandler ; 
#include "TH1F.h"

//_________________________________________________________________________ 
class AliAnalysisBtoJPSItoEle : public TNamed {
 public:
  //
  AliAnalysisBtoJPSItoEle();
  AliAnalysisBtoJPSItoEle(const AliAnalysisBtoJPSItoEle& source);
  AliAnalysisBtoJPSItoEle& operator=(const AliAnalysisBtoJPSItoEle& source);
  virtual ~AliAnalysisBtoJPSItoEle();
  
  Int_t DoMinimization(Double_t* x,Double_t* m, Int_t n);
  void ReadCandidates(TNtuple* nt, Double_t* &x, Double_t* &m, Int_t& n); // primary JPSI + secondary JPSI + bkg couples

  void SetPtBin(Int_t BinNum) { fPtBin = BinNum ; }
  void SetCsiMC(TH1F* MCtemplate);
  void SetResolutionConstants(Int_t BinNum);
  void CloneMCtemplate(const TH1F* MCtemplate) {fMCtemplate = (TH1F*)MCtemplate->Clone("fMCtemplate");}
  Double_t* GetResolutionConstants(Double_t* resolutionConst);
  Int_t GetPtBin() const { return fPtBin ; }

 private:
  //
  AliBtoJPSItoEleCDFfitHandler* fFCNfunction; //! pointer to the interface class
  Int_t fPtBin;                               // number of pt bin in which the analysis is performes
  TH1F* fMCtemplate;			      //! template of the MC distribution for the x distribution of the secondary J/psi

  ClassDef(AliAnalysisBtoJPSItoEle,0); // AliAnalysisBtoJPSItoEle class
};

#endif
