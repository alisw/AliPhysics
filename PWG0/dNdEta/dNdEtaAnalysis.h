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
  enum { kVertexBinning = 1+4 }; // the first is for the whole vertex range, the others divide the vertex range

  dNdEtaAnalysis(Char_t* name, Char_t* title);

  void FillTrack(Float_t vtx, Float_t eta);
  void FillEvent(Float_t vtx);

  void Finish(dNdEtaCorrection* correction);

  void DrawHistograms();
  void LoadHistograms();
  void SaveHistograms();

  virtual Long64_t Merge(TCollection* list);

  TH2F* GetEtaVsVtxHistogram() { return fEtaVsVtx; }
  TH2F* GetEtaVsVtxUncorrectedHistogram() { return fEtaVsVtxUncorrected; }
  TH1D* GetVtxHistogram() { return fVtx; }
  TH1D* GetdNdEtaHistogram(Int_t i = 0) { return fdNdEta[i]; }

protected:
  TH2F* fEtaVsVtx;              // histogram Eta vs vtx (track count)
  TH2F* fEtaVsVtxUncorrected;   // uncorrected histograms Eta vs vtx (track count)
  TH1D* fVtx;                   // vtx histogram (event count)
  TH1D* fdNdEta[kVertexBinning];// dndeta results for different vertex bins (0 = full range)

  ClassDef(dNdEtaAnalysis, 0)
};

#endif
