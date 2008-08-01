#ifndef ALIMUONSURVEYUTIL_H
#define ALIMUONSURVEYUTIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup geometry 
/// \class AliMUONSurveyUtil
/// \brief Utility class for survey of muon spectrometer
//
// Authors: Javier Castillo

class AliMUONGeometryTransformer;
class TGeoCombiTrans;
class TClonesArray;

class AliMUONSurveyUtil:public TObject
{

 public:
  /// Destructor
  virtual ~AliMUONSurveyUtil();

  static AliMUONSurveyUtil *Instance();

  static Bool_t MatrixToAngles(const Double_t *rot, Double_t *angles);
  static void AnglesToMatrix(const Double_t *angles, Double_t *rot);

  Double_t xpCenter(Double_t *x, Double_t *par);
  Double_t xnCenter(Double_t *x, Double_t *par);
  Double_t ypCenter(Double_t *x, Double_t *par);
  Double_t ynCenter(Double_t *x, Double_t *par);
  Double_t zpCenter(Double_t *x, Double_t *par);
  Double_t znCenter(Double_t *x, Double_t *par);
  Double_t phixpp(Double_t *x, Double_t *par);
  Double_t phixpn(Double_t *x, Double_t *par);
  Double_t phixnp(Double_t *x, Double_t *par);
  Double_t phixnn(Double_t *x, Double_t *par);
  Double_t phiypp(Double_t *x, Double_t *par);
  Double_t phiypn(Double_t *x, Double_t *par);
  Double_t phiynp(Double_t *x, Double_t *par);
  Double_t phiynn(Double_t *x, Double_t *par);
  
  static AliMUONGeometryTransformer *ReAlign(const AliMUONGeometryTransformer * transformer, 
					     int rMod, int rNDetElems, int rDetElemToDetElemId[], TGeoCombiTrans deltaDetElemTransf[], Bool_t verbose);

  static  void SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t chId, Double_t chResX, Double_t chResY, Double_t deResX, Double_t deResY);
  
 protected:   
  /// Default constructor
  AliMUONSurveyUtil() : TObject() {}

 private:
  static int fgNDetElemCh[10];
  static AliMUONSurveyUtil *fgInstance;

ClassDef(AliMUONSurveyUtil, 0) //Class for alignment of muon spectrometer
};

#endif
