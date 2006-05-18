#ifndef DNDETANALYSIS_H
#define DNDETANALYSIS_H

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
// - implement destructor

#include <TObject.h>
#include <TString.h>

class TH2F;
class TH1D;

class dNdEtaAnalysis : public TObject
{
public:
  dNdEtaAnalysis(Char_t* name="dndeta_correction");

  void FillTrack(Float_t vtx, Float_t eta, Float_t weight);
  void FillEvent(Float_t vtx);

  void Finish();

  void DrawHistograms();
  void SaveHistograms();

  TH2F* GetEtaVsVtxHistogram() { return hEtaVsVtx; }
  TH2F* GetEtaVsVtxUncorrectedHistogram() { return hEtaVsVtxUncorrected; }
  TH1D* GetVtxHistogram() { return hVtx; }
  TH1D* GetdNdEtaHistogram() { return hdNdEta; }

  void SetEtaVsVtxHistogram(TH2F* aHist) { hEtaVsVtx = aHist; }
  void SetEtaVsVtxUncorrectedHistogram(TH2F* aHist) { hEtaVsVtxUncorrected = aHist; }
  void SetVtxHistogram(TH1D* aHist) { hVtx = aHist; }
  void SetdNdEtaHistogram(TH1D* aHist) { hdNdEta = aHist; }

protected:
  TString  fName;

  TH2F* hEtaVsVtx;
  TH2F* hEtaVsVtxUncorrected;
  TH1D* hVtx;
  TH1D* hdNdEta;

  ClassDef(dNdEtaAnalysis,0)
};

#endif
