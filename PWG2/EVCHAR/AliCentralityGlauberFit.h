#ifndef ALICENTRALITYGLAUBERFIT_H
#define ALICENTRALITYGLAUBERFIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*   Origin: Alberica Toia, CERN, Alberica.Toia@cern.ch                   */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class to determine centrality percentiles from 1D distributions          // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliCentralityGlauberFit : public TObject {

 public:
  
  AliCentralityGlauberFit();
  virtual ~AliCentralityGlauberFit();

  void SetOutputFile(TString outrootfilename);
  void AddHisto(TString name);
  void MakeFits(TString infilename);
  Float_t GetPercentileCrossSection();
  void SetGlauberParam(Int_t fNmu, Float_t fmulow, Float_t fmuhigh, Int_t fNk, Float_t fklow, Float_t fkhigh, Int_t fNalpha, Float_t falphalow, Float_t falphahigh); 

 private:
  
  Int_t   Nmu;
  Float_t mulow;
  Float_t muhigh;
  Int_t   Nk;
  Float_t klow;
  Float_t khigh;
  Int_t   Nalpha;
  Float_t alphalow;
  Float_t alphahigh;

  TNtuple *glau_ntuple;
  Float_t npart;
  Float_t ncoll;
  Float_t b;

  Float_t effi;
  TH1D *heffi;

  TFile *inrootfile;
  TString outrootfilename;
  vector<TString> histnames;
  Float_t percentXsec;

  TH1D * NormalizeHisto(TString hdistributionName);
  TH1D * GlauberHisto(Float_t mu, Float_t k, Float_t eff, Float_t alpha, TH1D *hDATA, Bool_t save=kFALSE); 
  Float_t CalculateChi2(TH1D *hDATA, TH1D *thistGlau, Float_t eff);
  TH1D * GetTriggerEfficiencyFunction(TH1D *hist1, TH1D *hist2);
  Float_t GetTriggerEfficiencyIntegral(TH1D *hist1, TH1D *hist2); 

  Double_t NBD(Int_t n, Double_t mu, Double_t k);
  TH1D *NBDhist(Double_t mu, Double_t k);

  void  SaveHisto(TH1D *hist1,TH1D *hist2,TH1D *heffi, TFile *outrootfile);
  ClassDef(AliCentralityGlauberFit, 1)  
};
#endif


