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

#include "TObject.h"

// forward decl
class TF1;

class AliCentralityByFunction : public TObject {

 public:
  
  AliCentralityByFunction();
  virtual ~AliCentralityByFunction();

  void SetPercentileFile(TString outrootfilename);
  void SetPercentileCrossSection(Float_t percentXsec);
  //  void SetFitFunction(TString distribution, TString func);
  void SetFitFunction(TString distribution, TString func, Double_t xmin, Double_t xmax);
  void AddHisto(TString name);
  void MakePercentiles(TString infilename);

 private:

  TFile *inrootfile;
  TString outrootfilename;
  TFile *outrootfile;

  vector<TString> histnames;
  Float_t percentXsec;
  map<TString, TString>fitfunc;    // mapping from distribution to fit function name
  map<TString, TF1 *>fitter;  // mapping from fit function name to corresponding TF1
    
  TH1D *FitHisto(TString hdistributionName);
  TH1D * MakePercentHisto(TH2D *hist);

  ClassDef(AliCentralityByFunction, 1)  
};
#endif


