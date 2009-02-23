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
#include "AlidNdEtaCorrection.h"
#include "AliPWG0Helper.h"
#include <TString.h>

class TH1F;
class TCollection;

class AliCorrection;

class dNdEtaAnalysis : public TNamed
{
public:
  enum { kVertexBinning = 1+2 }; // the first is for the whole vertex range, the others divide the vertex range

  dNdEtaAnalysis();
  dNdEtaAnalysis(const Char_t* name, const Char_t* title, AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kSPD);
  virtual ~dNdEtaAnalysis();

  dNdEtaAnalysis(const dNdEtaAnalysis &c);
  dNdEtaAnalysis &operator=(const dNdEtaAnalysis &c);
  virtual void Copy(TObject &c) const;

  void FillTrack(Float_t vtx, Float_t eta, Float_t pt);
  void FillEvent(Float_t vtx, Float_t n);
  void FillTriggeredEvent(Float_t n);

  void Finish(AlidNdEtaCorrection* correction, Float_t ptCut, AlidNdEtaCorrection::CorrectionType correctionType, const char* tag = "");

  void DrawHistograms(Bool_t simple = kFALSE);
  void LoadHistograms(const Char_t* dir = 0);
  void SaveHistograms();

  virtual Long64_t Merge(TCollection* list);

  AliCorrection* GetData() const { return fData; }

  TH1F* GetPtHistogram() const { return fPtDist; }

  TH1F* GetdNdEtaHistogram(Int_t i = 0) const { return fdNdEta[i]; }
  TH1F* GetdNdEtaPtCutOffCorrectedHistogram(Int_t i = 0) const { return fdNdEtaPtCutOffCorrected[i]; }

  void SetAnalysisMode(AliPWG0Helper::AnalysisMode mode) { fAnalysisMode = mode; }
  AliPWG0Helper::AnalysisMode GetAnalysisMode() { return fAnalysisMode; }

protected:
  Float_t GetVtxMin(Float_t eta);

  AliCorrection* fData;     // we store the data in an AliCorrection
  TH1F* fMult;              // (raw) multiplicity distribution of triggered events

  TH1F* fPtDist; // pt distribution

  TH1F* fdNdEta[kVertexBinning];                   // dndeta results for different vertex bins (0 = full range)
  TH1F* fdNdEtaPtCutOffCorrected[kVertexBinning];  // dndeta results for different vertex bins (0 = full range), pt cut off corrected

  AliPWG0Helper::AnalysisMode fAnalysisMode;       // detector (combination) used for the analysis
  TString fTag;                                   // tag saved that describes the applied correction

  ClassDef(dNdEtaAnalysis, 2)
};

#endif
