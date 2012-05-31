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

  Double_t XpCenter(const Double_t *x, const Double_t *par) const;
  Double_t XnCenter(const Double_t *x, const Double_t *par) const;
  Double_t YpCenter(const Double_t *x, const Double_t *par) const;
  Double_t YnCenter(const Double_t *x, const Double_t *par) const;
  Double_t ZpCenter(const Double_t *x, const Double_t *par) const;
  Double_t ZnCenter(const Double_t *x, const Double_t *par) const;
  Double_t PhiXpp(const Double_t *x, const Double_t *par) const;
  Double_t PhiXpn(const Double_t *x, const Double_t *par) const;
  Double_t PhiXnp(const Double_t *x, const Double_t *par) const;
  Double_t PhiXnn(const Double_t *x, const Double_t *par) const;
  Double_t PhiYpp(const Double_t *x, const Double_t *par) const;
  Double_t PhiYpn(const Double_t *x, const Double_t *par) const;
  Double_t PhiYnp(const Double_t *x, const Double_t *par) const;
  Double_t PhiYnn(const Double_t *x, const Double_t *par) const;
  
  static AliMUONGeometryTransformer *ReAlign(const AliMUONGeometryTransformer * transformer, 
					     int rMod, int rNDetElems, int rDetElemToDetElemId[], TGeoCombiTrans deltaDetElemTransf[], Bool_t verbose);

  static  void SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t chId, Double_t chResX, Double_t chResY, Double_t deResX, Double_t deResY);
  
 protected:   
  /// Default constructor
  AliMUONSurveyUtil() : TObject() {}

 private:
  /// Not implemented
  AliMUONSurveyUtil(const AliMUONSurveyUtil& right);
  /// Not implemented
  AliMUONSurveyUtil&  operator = (const AliMUONSurveyUtil& right);


  static int fgNDetElemCh[10];  ///< Numbers of detection elements per chamber
  static AliMUONSurveyUtil *fgInstance;   ///< Singleton instance 

ClassDef(AliMUONSurveyUtil, 0) //Class for alignment of muon spectrometer
};

#endif
