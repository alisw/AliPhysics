#ifndef ALICORRECTION_H
#define ALICORRECTION_H

/* $Id$ */

// ------------------------------------------------------
//
// This class is used to store the correction for one effect. 
//
// Most effects have to be correction on track and event level, this class combines
// two correction matrices. One of the type AliCorrectionMatrix2D and one of
// the type AliCorrectionMatrix3D
//
// ------------------------------------------------------

#include <TNamed.h>
#include "AliPWG0Helper.h"

class AliCorrectionMatrix2D;
class AliCorrectionMatrix3D;

class AliCorrection : public TNamed
{
public:
  AliCorrection();
  AliCorrection(const Char_t* name, const Char_t* title, AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kTPC);
  AliCorrection(const AliCorrection& c);

  virtual ~AliCorrection();
  AliCorrection& operator=(const AliCorrection& corr);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);

  AliCorrectionMatrix2D* GetEventCorrection() const { return fEventCorr; }
  AliCorrectionMatrix3D* GetTrackCorrection() const { return fTrackCorr; }

  void SetEventCorrection(AliCorrectionMatrix2D* corr) { fEventCorr = corr; }
  void SetTrackCorrection(AliCorrectionMatrix3D* corr) { fTrackCorr = corr; }

  void Divide();
  void Multiply();
  void SetCorrectionToUnity();
  void Scale(Double_t factor);

  void Add(AliCorrection* aCorrectionToAdd, Float_t c=1);

  virtual Bool_t LoadHistograms(const Char_t* dir = 0);
  virtual void SaveHistograms();
  virtual void DrawHistograms(const Char_t* name = 0);
  virtual void DrawOverview(const char* canvasName = 0);

  virtual void ReduceInformation();

  virtual void Reset(Option_t* option = "");
  void PrintStats(Float_t zRange, Float_t etaRange, Float_t ptCut);
  void PrintInfo(Float_t ptCut);

protected:
  AliCorrectionMatrix2D* fEventCorr; // correction on event level
  AliCorrectionMatrix3D* fTrackCorr; // correction on track level

  ClassDef(AliCorrection,1)
};

#endif

