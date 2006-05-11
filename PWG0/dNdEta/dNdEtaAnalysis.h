// ------------------------------------------------------
//
// Class for dn/deta analysis
//
// ------------------------------------------------------
// 
// TODO:
// - more documentation
// - add debug statements
// - add more histograms
// - add functionality to set the bin sizes
// - figure out correct way to treat the errors
// - add functionality to make dn/deta for different mult classes?

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TFile
#include "TFile.h"
#endif
#ifndef ROOT_TH2
#include "TH2.h"
#endif
#ifndef ROOT_TMath
#include "TMath.h"
#endif

class dNdEtaAnalysis : public TObject 
{
protected:    

  TString  fName; 

  TH2F* hEtaVsVtx;
  TH1F* hVtx;
  TH1F* hdNdEta;

public:
  dNdEtaAnalysis(Char_t* name="dndeta_correction");

  void FillTrack(Float_t vtx, Float_t eta, Float_t weight);
  void FillEvent(Float_t vtx);
  
  void Finish();
  
  void SaveHistograms();
  
  ClassDef(dNdEtaAnalysis,0)
};


