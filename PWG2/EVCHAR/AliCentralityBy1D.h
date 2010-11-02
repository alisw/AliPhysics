#ifndef ALICENTRALITYBY1D_H
#define ALICENTRALITYBY1D_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*   Origin: Alberica Toia, CERN, Alberica.Toia@cern.ch                   */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class to determine centrality percentiles from 1D distributions          // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliCentralityBy1D : public TObject {

 public:
  
  AliCentralityBy1D();
  virtual ~AliCentralityBy1D();

  void SetPercentileFile(TString outrootfilename);
  void SetPercentileCrossSection(Float_t percentXsec);
  void AddHisto(TString name);
  void MakePercentiles(TString infilename);

 private:
  
  TFile *inrootfile;

  TString outrootfilename;
  vector<TString> histnames;
  Float_t percentXsec;

  TH1D * MakePercentHisto(TString hdistributionName);
  void  SaveHisto(TH1D *hist, TFile *outrootfile);

  ClassDef(AliCentralityBy1D, 1)  
};
#endif


