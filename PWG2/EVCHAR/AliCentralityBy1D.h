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

  void SetInputFile(TString filename)                 { fInrootfilename = filename;   }
  void SetOutputFile(TString filename)                { fOutrootfilename = filename;  }

  void SetPercentileCrossSection(Float_t xsec)        { fPercentXsec = xsec;  }
  void SetMultLowBound(Float_t mult )                 { fMultLowBound = mult; }

  void AddHisto(TString name) { fHistnames.push_back(name); }
  void MakePercentiles();

 private:
  TString fInrootfilename;          // input root file
  TString fOutrootfilename;         // output root file
  std::vector<TString> fHistnames;  // histogram names
  
  Float_t fPercentXsec;
  Float_t fMultLowBound;

  void  SaveHisto(TH1F *hist1, TH1F *hist2, TFile *outrootfile);

  ClassDef(AliCentralityBy1D, 1)  
};
#endif


