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

class AliCorrectionMatrix2D;
class AliCorrectionMatrix3D;

class AliCorrection : public TNamed
{
public:
  AliCorrection();
  AliCorrection(const Char_t* name, const Char_t* title);
  AliCorrection(const AliCorrection& c);

  virtual ~AliCorrection();
  AliCorrection& operator=(const AliCorrection& corr);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);

  AliCorrectionMatrix2D* GetEventCorrection() { return fEventCorr; }
  AliCorrectionMatrix3D* GetTrackCorrection() { return fTrackCorr; }

  void SetEventCorrection(AliCorrectionMatrix2D* corr) { fEventCorr = corr; }
  void SetTrackCorrection(AliCorrectionMatrix3D* corr) { fTrackCorr = corr; }

  void Divide();
  void Multiply();
  void SetCorrectionToUnity();

  void Add(AliCorrection* aCorrectionToAdd, Float_t c=1);

  virtual Bool_t LoadHistograms(const Char_t* dir = 0);
  virtual void SaveHistograms();
  virtual void DrawHistograms(const Char_t* name = 0);

  virtual void ReduceInformation();

  virtual void Reset(Option_t* option = "");

protected:
  AliCorrectionMatrix2D* fEventCorr; // correction on event level
  AliCorrectionMatrix3D* fTrackCorr; // correction on track level

  ClassDef(AliCorrection,1)
};

#endif

