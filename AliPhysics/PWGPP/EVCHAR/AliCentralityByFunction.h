#ifndef ALICENTRALITYBYFUNCTION_H
#define ALICENTRALITYBYFUNCTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*   Origin: Alberica Toia, CERN, Alberica.Toia@cern.ch                   */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class to determine centrality percentiles from 2D distributions          // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>
#include <TString.h>

// forward decl
class TF1;
class TH1D;
class TH2D;

class AliCentralityByFunction : public TObject {

 public:
  
  AliCentralityByFunction();
  virtual ~AliCentralityByFunction() {}

  void SetPercentileFile(TString filename)            { foutrootfilename = filename; }
  void SetPercentileCrossSection(Float_t xsec)        { fpercentXsec     = xsec;     }
  void SetFitFunction(TString distribution, TString func, Double_t xmin, Double_t xmax);
  void AddHisto(TString name)                         { fhistnames.push_back(name); }
  void MakePercentiles(TString infilename);

 private:
  AliCentralityByFunction(const AliCentralityByFunction&);
  AliCentralityByFunction& operator=(const AliCentralityByFunction&);

  TFile   *finrootfile;               // input root file
  TString  foutrootfilename;          // output root file name
  TFile   *foutrootfile;              // output root file
  std::vector<TString>    fhistnames; // hist names
  Float_t               fpercentXsec; // percentile cross section
  std::map<TString, TString> fitfunc; // mapping from distribution to fit function name
  std::map<TString, TF1 *>   fitter;  // mapping from fit function name to corresponding TF1
    
  TH1D *FitHisto(TString hdistributionName);
  TH1D *MakePercentHisto(TH2D *hist);

  ClassDef(AliCentralityByFunction, 1)  
};
#endif
