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
  Double_t GetPercentileCrossSection() { return fPercentXsec;         }
  TH1F    *GetTempHist()               { return fTempHist;            } 
  void     MakeFits();
  void     MakeFitsMinuitNBD(Double_t alpha=-1, Double_t mu=-1, Double_t k=-1, Double_t eff=-1);
  void     SetAncestorMode(Int_t f)  { fAncestor=f; }
  void     SetFastFitMode(Int_t f)  { fFastFit=f; }
  void     SetGlauberParam(Int_t Nmu, Double_t mulow, Double_t muhigh, Int_t Nk, Double_t klow, 
                          Double_t khigh, Int_t Nalpha, Double_t alphalow, Double_t alphahigh, 
                          Int_t Neff, Double_t efflow, Double_t effhigh); 
  void     SetInputFile(TString filename)         { fInrootfilename = filename;   }
  void     SetInputNtuple(TString ntuplename)     { fInntuplename   = ntuplename; }
  void     SetAncFile(TString name)               { fAncfilename    = name;       }
  void     SetNevents(Int_t f) { fNevents=f; }
  void     SetNtrials(Int_t f) { fNtrials=f; }
  void     SetOutputFile(TString filename)        { fOutrootfilename = filename;  }
  void     SetOutputNtuple(TString ntuplename)    { fOutntuplename   = ntuplename;}
  void     SetRangeToFit(Double_t multmin, Double_t multmax);
  void     SetRebin(Int_t f)  { fRebinFactor=f; }
  void     UseAverage(Bool_t f)  { fUseAverage=f; }
  void     UseChi2(Bool_t f)  { fUseChi2=f; }

 private:
  Int_t   fNmu;           // mu points
  Double_t fMulow;        // mu low 
  Double_t fMuhigh;       // mu high
  Int_t   fNk;            // k points
  Double_t fKlow;         // k low
  Double_t fKhigh;        // k high
  Int_t   fNalpha;        // alpha points
  Double_t fAlphalow;     // alpha low
  Double_t fAlphahigh;    // alpha high
  Int_t   fNeff;          // eff points
  Double_t fEfflow;       // eff low
  Double_t fEffhigh;      // eff high
  Int_t   fRebinFactor;   // rebin factor
  Double_t fMultmin;      // mult min
  Double_t fMultmax;      // mult max
  TNtuple *fGlauntuple;   // glauber ntuple
  Float_t fNpart;         // number participants
  Float_t fNcoll;         // number collisions
  Float_t fB;             // impact parameter
  Float_t fTaa;           // taa
  Float_t fEffi;          // efficiency
  TH1F   *fhEffi;         // efficiency histogram
  TH1F  *fTempHist;       // Temporary pointer to data histo, to be used in minimization 
  TH1F  *fGlauberHist;    // Glauber histogram
  Int_t fFastFit;         // If 1 or 2 use cruder approximation to compute curve faster 1:NBD, 2:Gauss
  Int_t fAncestor;        // If 1 use Npart**alpha, if 2 use alpha*Npart + (1-alpha)*Ncoll
  TF1 * fNBD;             // NBD function
  Bool_t fUseChi2;        // If true, use chi2
  Bool_t fUseAverage;     // If true, use average
  TH1F *fhAncestor;       // histo for the ancestor distribution
  Int_t fNevents;         // Number of events to use in the glauber ntuple
  Int_t fNtrials;         // Number of trials to use for the average

  TString fInrootfilename;          // input root file
  TString fInntuplename;            // input Gauber ntuple
  TString fOutrootfilename;         // output root file
  TString fOutntuplename;           // output Glauber ntuple
  TString fAncfilename;             // ancestor file name
  std::vector<TString> fHistnames;  // histogram names
  Double_t fPercentXsec;            // cross section 

  Double_t  CalculateChi2(TH1F *hDATA, TH1F *thistGlau, Double_t eff);
  TH1F     *GetTriggerEfficiencyFunction(TH1F *hist1, TH1F *hist2);
  Double_t  GetTriggerEfficiencyIntegral(TH1F *hist1, TH1F *hist2); 
  TH1F     *GlauberHisto(Double_t mu, Double_t k, Double_t eff, Double_t alpha, TH1F *hDATA, Bool_t save=kFALSE); 
  TH1F     *MakeAncestor(Double_t alpha);
  Double_t  NBD(Int_t n, Double_t mu, Double_t k);
  TH1F     *NBDhist(Double_t mu, Double_t k);
  TH1F     *NormalizeHisto(TString hdistributionName);
  void      SaveHisto(TH1F *hist1,TH1F *hist2,TH1F *heffi, TFile *outrootfile);

  static void     MinuitFcnNBD(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  static Double_t NBDFunc(Double_t *p, Double_t * x);

  AliCentralityGlauberFit(const AliCentralityGlauberFit&);
  AliCentralityGlauberFit &operator=(const AliCentralityGlauberFit&);

  ClassDef(AliCentralityGlauberFit, 1)  
};
#endif
