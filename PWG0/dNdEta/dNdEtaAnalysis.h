/* $Id$ */

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

#include <TNamed.h>

class TH2F;
class TH1D;
class TCollection;
class dNdEtaCorrection;

class dNdEtaAnalysis : public TNamed
{
public:
  enum { kVertexBinning = 1+6 }; // the first is for the whole vertex range, the others divide the vertex range

  dNdEtaAnalysis(Char_t* name, Char_t* title);

  void FillTrack(Float_t vtx, Float_t eta, Float_t c);
  void FillEvent(Float_t vtx);

  void Finish(dNdEtaCorrection* correction);

  void DrawHistograms();
  void SaveHistograms();

  virtual Long64_t Merge(TCollection* list);

  TH2F* GetEtaVsVtxHistogram() { return hEtaVsVtx; }
  TH2F* GetEtaVsVtxUncorrectedHistogram() { return hEtaVsVtxUncorrected; }
  TH1D* GetVtxHistogram() { return hVtx; }
  TH1D* GetdNdEtaHistogram(Int_t i = 0) { return hdNdEta[i]; }

protected:
  TH2F* hEtaVsVtx;
  TH2F* hEtaVsVtxCheck;
  TH2F* hEtaVsVtxUncorrected;
  TH1D* hVtx;
  TH1D* hdNdEta[kVertexBinning];

  ClassDef(dNdEtaAnalysis,0)
};

#endif
