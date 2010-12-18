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

#include <vector>
#include <TObject.h>
#include <TString.h>

class TNtuple;
class TFile;

//using namespace std;


class AliCentralityGlauberFit : public TObject {

 public:
  
  AliCentralityGlauberFit(const char * filename="/home/alberica/GlauberNtuple/GlauberMC_PbPb_ntuple_sigma64_mind4_r662_a546.root");
  virtual ~AliCentralityGlauberFit() {}

  void    SetOutputFile(TString filename)        { foutrootfilename = filename; }
  void    AddHisto(TString name)                 { fhistnames.push_back(name);  }
  void    MakeFits(TString infilename);
  void    MakeFitsMinuitNBD(TString infilename);
  Float_t GetPercentileCrossSection() { return fPercentXsec; }
  void    SetGlauberParam(Int_t Nmu, Float_t mulow, Float_t muhigh, Int_t Nk, Float_t klow, 
                          Float_t khigh, Int_t Nalpha, Float_t alphalow, Float_t alphahigh, 
                          Int_t Neff, Float_t efflow, Float_t effhigh); 
  void    SetRangeToFit(Float_t multmin, Float_t multmax);
  void    SetRebin(Int_t f)  { fRebinFactor=f; }
  TH1D   *GetTempHist()    { return fTempHist;    } 
  TH1D   *GetGlauberHist() { return fGlauberHist; } 
  static void MinuitFcnNBD(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

 private:
  Int_t   fNmu;          // mu points
  Float_t fMulow;        // mu low 
  Float_t fMuhigh;       // mu high
  Int_t   fNk;           // k points
  Float_t fKlow;         // k low
  Float_t fKhigh;        // k high
  Int_t   fNalpha;       // alpha points
  Float_t fAlphalow;     // alpha low
  Float_t fAlphahigh;    // alpha high
  Int_t   fNeff;         // eff points
  Float_t fEfflow;       // eff low
  Float_t fEffhigh;      // eff high
  Int_t   fRebinFactor;  // rebin factor
  Float_t fMultmin;      // mult min
  Float_t fMultmax;      // mult max
  TNtuple *fGlauntuple;  // glauber ntuple
  Float_t fNpart;        // number participants
  Float_t fNcoll;        // number collisions
  Float_t fB;            // impact parameter
  Float_t fTaa;          // taa
  Float_t fEffi;         // efficiency
  TH1D   *fhEffi;        // efficiency histogram
  TH1D  *fTempHist;      // Temporary pointer to data histo, to be used in minimization 
  TH1D  *fGlauberHist;   // Glauber histogram
  Int_t fFastFit;        // If true, use cruder approximation to compute curve faster
  TF1 * fNBD;            // NBD function
  Bool_t fUseChi2;       // If true, use chi2

  TFile *finrootfile;              // input root file
  TString foutrootfilename;        // output root file
  std::vector<TString> fhistnames; // histogram names
  Float_t fPercentXsec;            // cross section 

  TH1D    *NormalizeHisto(TString hdistributionName);
  TH1D    *GlauberHisto(Float_t mu, Float_t k, Float_t eff, Float_t alpha, TH1D *hDATA, Bool_t save=kFALSE); 
  Float_t  CalculateChi2(TH1D *hDATA, TH1D *thistGlau, Float_t eff);
  TH1D    *GetTriggerEfficiencyFunction(TH1D *hist1, TH1D *hist2);
  Float_t  GetTriggerEfficiencyIntegral(TH1D *hist1, TH1D *hist2); 
  static Double_t NBDFunc(Double_t *p, Double_t * x);
  Double_t NBD(Int_t n, Double_t mu, Double_t k);
  TH1D    *NBDhist(Double_t mu, Double_t k);
  void     SaveHisto(TH1D *hist1,TH1D *hist2,TH1D *heffi, TFile *outrootfile);

  AliCentralityGlauberFit(const AliCentralityGlauberFit&);
  AliCentralityGlauberFit &operator=(const AliCentralityGlauberFit&);

  ClassDef(AliCentralityGlauberFit, 1)  
};
#endif
