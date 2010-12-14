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
#include <vector>
#include "TString.h"

using namespace std;


class AliCentralityGlauberFit : public TObject {

 public:
  
  AliCentralityGlauberFit(const char * filename="../centrality/files/GlauberMC_PbPb_ntuple_sigma64_mind4_r662_a546.root");
  virtual ~AliCentralityGlauberFit();

  void SetOutputFile(TString outrootfilename);
  void AddHisto(TString name);
  void MakeFits(TString infilename);
  void MakeFitsMinuitNBD(TString infilename);
  Float_t GetPercentileCrossSection();
  void SetGlauberParam(Int_t fNmu, Float_t fmulow, Float_t fmuhigh, Int_t fNk, Float_t fklow, Float_t fkhigh, Int_t fNalpha, Float_t falphalow, Float_t falphahigh, Int_t fNeff, Float_t fefflow, Float_t feffhigh); 
  void SetRangeToFit(Float_t fmultmin, Float_t fmultmax);
  void SetRebin(Int_t frebinFactor);
  static void MinuitFcnNBD(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  TH1D * GetTempHist() { return fTempHist;} 
  TH1D * GetGlauberHist() { return fGlauberHist;} 


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
  Int_t   Neff;
  Float_t efflow;
  Float_t effhigh;

  Int_t rebinFactor;
  Float_t multmin;
  Float_t multmax;

  TNtuple *glau_ntuple;
  Float_t Npart;
  Float_t Ncoll;
  Float_t B;
  Float_t tAA;

  Float_t effi;
  TH1D *heffi;
  
  TH1D * fTempHist;    // Temporary pointer to data histo, to be used in minimization 
  TH1D * fGlauberHist; // Glauber histogram
  Int_t fFastFit; // If true, use cruder approximation to compute curve faster
  TF1 * fNBD; // NBD function
  Bool_t fUseChi2;

  TFile *inrootfile;
  TString outrootfilename;
  vector<TString> histnames;
  Float_t percentXsec;

  TH1D * NormalizeHisto(TString hdistributionName);
  TH1D * GlauberHisto(Float_t mu, Float_t k, Float_t eff, Float_t alpha, TH1D *hDATA, Bool_t save=kFALSE); 
  Float_t CalculateChi2(TH1D *hDATA, TH1D *thistGlau, Float_t eff);
  TH1D * GetTriggerEfficiencyFunction(TH1D *hist1, TH1D *hist2);
  Float_t GetTriggerEfficiencyIntegral(TH1D *hist1, TH1D *hist2); 

  static Double_t NBDFunc(Double_t *p, Double_t * x);

  Double_t NBD(Int_t n, Double_t mu, Double_t k);
  TH1D *NBDhist(Double_t mu, Double_t k);

  void  SaveHisto(TH1D *hist1,TH1D *hist2,TH1D *heffi, TFile *outrootfile);
  ClassDef(AliCentralityGlauberFit, 1)  
};
#endif


