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

class AliCentralityGlauberFit : public TObject {

 public:
  
  AliCentralityGlauberFit(const char * filename);
  virtual ~AliCentralityGlauberFit() {}

  void     AddHisto(TString name)      { fHistnames.push_back(name);  }
  TH1F    *GetGlauberHist()            { return fGlauberHist;         } 
  void     MakeFits();
  void     SetGlauberParam(Int_t Nk, Double_t klow, Double_t khigh, Int_t Nalpha, Double_t alphalow, Double_t alphahigh, 
  		Int_t Nsigma, Double_t sigmalow, Double_t sigmahigh, Int_t Nbog, Double_t boglow, Double_t boghigh, 
		Int_t Ngamma, Double_t gammalow, Double_t gammahigh, Int_t Nsigmap, Double_t sigmaplow, Double_t sigmaphigh); 
  void     SetNParam(double a, double b, double c, double d); 
  void     SetInputFile(TString filename)         { fInrootfilename = filename;   }
  void     SetInputNtuple(TString ntuplename)     { fInntuplename   = ntuplename; }
  void     SetNevents(Int_t f) { fNevents=f; }
  void     SetOutputFile(TString filename)        { fOutrootfilename = filename;  }
  void     SetOutputNtuple(TString ntuplename)    { fOutntuplename   = ntuplename;}
  void     SetRangeToFit(Double_t multmin, Double_t multmax);
  void     SetRangeToScale(Double_t scalemin);
  void     SetRebin(Int_t f)  { fRebinFactor=f; }
  void     UseChi2(Bool_t f)  { fUseChi2=f; }
  void     SetSaturation(Bool_t saturation) {fApplySaturation = saturation;}
  void     SetSaturationParams(Int_t ngray=15, Int_t nblack=28) {fnGraySaturation=ngray; fnBlackSaturation=nblack;}
  void     SetIsZN() {fIsZN=kTRUE;}
  void     SetIsZP() {fIsZN=kFALSE;}
  void     SetNTrials(int ntrials) {fNTrials=ntrials;}
  void     SetChi2Min(double chi2min) {fChi2Min=chi2min;}
  
 private:
  Int_t    fNk;           // k points
  Double_t fKlow;         // k low
  Double_t fKhigh;        // k high
  Int_t    fNalpha;       // alpha points
  Double_t fAlphalow;     // alpha low
  Double_t fAlphahigh;    // alpha high
  Int_t    fNsigma;       // sigma points
  Double_t fSigmalow;     // sigma low
  Double_t fSigmahigh;    // sigma high
  Int_t    fNbog;         // bog points
  Double_t fBoglow;       // bog low
  Double_t fBoghigh;      // bog high
  Int_t    fNgamma;       // gamma points
  Double_t fgammalow;     // gamma low
  Double_t fgammahigh;    // gamma high
  Int_t    fNsigmap;      // sigmap points
  Double_t fsigmaplow;    // sigmap low
  Double_t fsigmaphigh;   // sigmap high
  Int_t    fRebinFactor;  // rebin factor
  Double_t fScalemin;     // mult min where to scale
  Double_t fMultmin;      // mult min
  Double_t fMultmax;      // mult max
  TNtuple *fGlauntuple;   // glauber ntuple
  Float_t  fNpart;        // number participants
  Float_t  fNcoll;        // number collisions
  Float_t  fB;            // impact parameter
  Float_t  fTaa;          // taa
  TH1F    *fTempHist;     // Temporary pointer to data histo, to be used in minimization 
  TH1F    *fGlauberHist;     // Glauber histogram
  Bool_t   fUseChi2;         // If true, use chi2
  Int_t    fNevents;         // Number of events to use in the glauber ntuple
  Bool_t   fApplySaturation; 
  Int_t    fnGraySaturation; 
  Int_t    fnBlackSaturation;
  Bool_t   fIsZN;
  Float_t  fSDprobability[40];
  Double_t fParamnSlow[4];
  Int_t    fNTrials;
  Double_t fChi2Min;

  TString fInrootfilename;          // input root file
  TString fInntuplename;            // input Gauber ntuple
  TString fOutrootfilename;         // output root file
  TString fOutntuplename;           // output Glauber ntuple
  TString fAncfilename;             // ancestor file name
  std::vector<TString> fHistnames;  // histogram names

  Double_t CalculateChi2(TH1F *hDATA, TH1F *thistGlau);
  TH1F    *GlauberHisto(Double_t k, Double_t alpha, Double_t sigma, Double_t bog, Double_t gamma, Double_t sigmap,
  		TH1F *hDATA, Bool_t save=kFALSE); 
  void     SaveHisto(TH1F *hist1,TH1F *hist2, TFile *outrootfile);
  void     MakeSlowNucleons(Int_t fNcoll, Double_t &nbn, Double_t &ngn,Double_t &nbp, Double_t &ngp);
  void     MakeSlowNucleons2(Int_t fNcoll, Double_t bog, Double_t gamma, 
  		Double_t &nbn, Double_t &ngn,Double_t &nbp, Double_t &ngp, Double_t &lcp);
  void     MakeSlowNucleons2s(Float_t fNcoll, Double_t bog, Double_t gamma, Double_t sigmap,
  		Double_t &nbn, Double_t &ngn,Double_t &nbp, Double_t &ngp, Double_t &lcp);
  Double_t ConvertToEnergy(Double_t T);
  Double_t Maxwell(Double_t m, Double_t p, Double_t T);

  AliCentralityGlauberFit(const AliCentralityGlauberFit&);
  AliCentralityGlauberFit &operator=(const AliCentralityGlauberFit&);

  ClassDef(AliCentralityGlauberFit, 4)  
};
#endif
