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

#include <TNamed.h>

class TH3F;
class TH1D;
class TCollection;
class AlidNdEtaCorrection;

class dNdEtaAnalysis : public TNamed
{
public:
  enum { kVertexBinning = 1+3 }; // the first is for the whole vertex range, the others divide the vertex range

  dNdEtaAnalysis();
  dNdEtaAnalysis(Char_t* name, Char_t* title);
  virtual ~dNdEtaAnalysis();

  dNdEtaAnalysis(const dNdEtaAnalysis &c);
  dNdEtaAnalysis &operator=(const dNdEtaAnalysis &c);
  virtual void Copy(TObject &c) const;

  void FillTrack(Float_t vtx, Float_t eta, Float_t pt, Float_t weight);
  void FillEvent(Float_t vtx, Float_t weight);

  void Finish(AlidNdEtaCorrection* correction, Float_t ptCut);

  void DrawHistograms();
  void LoadHistograms();
  void SaveHistograms();

  virtual Long64_t Merge(TCollection* list);

  TH3F* GetHistogram() { return fData; }
  TH3F* GetUncorrectedHistogram() { return fDataUncorrected; }
  TH1D* GetVtxHistogram() { return fVtx; }
  TH1D* GetdNdEtaHistogram(Int_t i = 0) { return fdNdEta[i]; }

protected:
  TH3F* fData;              // histogram Eta vs vtx (track count)
  TH3F* fDataUncorrected;   // uncorrected histograms Eta vs vtx (track count)

  TH1D* fVtx;                   // vtx histogram (event count)

  TH1D* fPtDist; // pt distribution

  TH1D* fdNdEta[kVertexBinning]; // dndeta results for different vertex bins (0 = full range)
  TH1D* fdNdEtaPtCutOffCorrected[kVertexBinning];  // dndeta results for different vertex bins (0 = full range), pt cut off corrected

  ClassDef(dNdEtaAnalysis, 1)
};

#endif
