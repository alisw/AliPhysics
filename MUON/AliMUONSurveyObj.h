#ifndef ALIMUONSURVEYOBJ_H
#define ALIMUONSURVEYOBJ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup geometry
/// \class AliMUONSurveyObj
/// \brief Base class for survey of muon spectrometer
//
// Author: Javier Castillo

#include <TObject.h>

class TObjArray;
class TGeoCombiTrans;
class TVector3;
class TH2;
class TF2;
class TFitter;
class TArrayD;

class AliSurveyPoint;

class AliMUONSurveyObj:public TObject
{

 public:
  AliMUONSurveyObj();
  virtual ~AliMUONSurveyObj();
 
  virtual Int_t AddStickerTargets(TObjArray *pArray, TString stBaseName, Int_t lTargetMax = 9);
  virtual Int_t AddGButtonTargets(TObjArray *pArray, TString btBaseName, Int_t lTargetMax = 9);
  virtual Int_t AddLButtonTargets(TObjArray *pArray, TString btBaseName, Int_t lTargetMax = 9);

  virtual Int_t AddStickerTargets(TString stBaseName, Int_t lTargetMax = 9) = 0;
  virtual Int_t AddGButtonTargets(TString btBaseName, Int_t lTargetMax = 9) = 0;

  void AddStickerTarget(AliSurveyPoint *stPoint);
  void AddGButtonTarget(AliSurveyPoint *btPoint);
  void AddLButtonTarget(AliSurveyPoint *btPoint);
  void AddLButtonTarget(TVector3 *btVector);

  Int_t GetNStickerTargets();
  AliSurveyPoint *GetStickerTarget(Int_t stIndex);
  Int_t GetNGButtonTargets();
  AliSurveyPoint *GetGButtonTarget(Int_t btIndex);
  Int_t GetNLButtonTargets();
  AliSurveyPoint *GetLButtonTarget(Int_t btIndex);
  /// Set transformation of geoemtrical element
  void SetBaseTransformation(TGeoCombiTrans *baseTrf, Bool_t ownerBaseTrf = kFALSE) {
    fBaseTrf = baseTrf;
    fOwnerBaseTrf=ownerBaseTrf;
  }
  /// Set local transformation of geometrical element
  virtual void SetLocalTransformation(TGeoCombiTrans *localTrf, Bool_t ownerLocalTrf = kFALSE) {
    fLocalTrf=localTrf;
    fOwnerLocalTrf=ownerLocalTrf;
  }

  TGeoCombiTrans* GetLocalTrf() const {return fLocalTrf;} 
  TGeoCombiTrans* GetBaseTrf() const {return fBaseTrf;}
  TGeoCombiTrans* GetAlignTrf()const {return fAlignTrf;}

  void SetUseCM(Bool_t bUseCM = kTRUE) {fUseCM = bUseCM;}
  Bool_t GetUseCM() const {return fUseCM;}

  void SetPlane(TString pName, Double_t xMin=-2000., Double_t xMax=+2000., Double_t yMin=-2000., Double_t yMax=2000.);
  void SetPlaneParameters(Double_t p0, Double_t p1, Double_t p2);

  void DrawSTargets();
  Double_t FitPlane();

  TF2* GetPlane() const {return fPlane;}

  TFitter* GetFitter() const {return fFitter;}

  Int_t SurveyToAlign(TGeoCombiTrans &quadTransf, Double_t *parErr, Double_t psi=0., Double_t tht=0., Double_t epsi=0., Double_t etht=0.);
  Int_t SurveyToAlign(Double_t psi=0., Double_t tht=0., Double_t epsi=0., Double_t etht=0.);
  Double_t SurveyChi2(Double_t *par);

  Double_t EvalFunction(const TF2 *lFunction, Int_t iP1, Int_t iP2, const Char_t *lCoord);

  void CalculateTranslation(TF2 *xFunc, TF2 *yFunc, TF2 *zFunc, Int_t iP1, Int_t iP2, Double_t *lCenTemp);
  //  TGeoCombiTrans *CalculateTransformation(TF2 *xFunc, TF2 *yFunc, TF2 *zFunc, TF2 *pFunc, Int_t iP1, Int_t iP2);

  Double_t CalculateGlobalDiff(TGeoCombiTrans &lTransf, Int_t nPoints, TArrayD &lDiff);

  Int_t CalculateBestTransf(Int_t iP1, Int_t iP2, Double_t *lXYZ, Double_t *lPTP);

  void CalculateMeanTransf(Double_t *lXYZ, Double_t *lPTP);

  void SetXMin(Double_t xMin) {fXMin = xMin;}
  void SetXMax(Double_t xMax) {fXMax = xMax;}
  void SetYMin(Double_t yMin) {fYMin = yMin;}
  void SetYMax(Double_t yMax) {fYMax = yMax;}
  void SetZMin(Double_t zMin) {fZMin = zMin;}
  void SetZMax(Double_t zMax) {fZMax = zMax;}

  virtual void PrintLocalTrf();
  virtual void PrintAlignTrf();

  void FillSTHistograms(TString baseNameC, TH2 *hSTc, TString baseNameA="", TH2 *hSTa = 0);

  Double_t GetAlignResX();
  Double_t GetAlignResY();

  AliSurveyPoint* ConvertPointUnits(AliSurveyPoint *stPoint, Float_t lFactor = 0.1);

 private:
  /// Not implemented
  AliMUONSurveyObj(const AliMUONSurveyObj& right);
  /// Not implemented
  AliMUONSurveyObj&  operator = (const AliMUONSurveyObj& right);

  Double_t EqPlane(const Double_t *x, const Double_t *par) const {
    return (-par[1]*x[0] +par[0]*x[1] -par[2]);  // then psi=ATan(par[0]) and tht=ATan(par[0])
    //    return (-par[0]*x[0] -par[1]*x[1] -par[2]); 
  }

  TObjArray *fSTargets;   ///< Array of AliSurveyPoint of Sticker Targets
  TObjArray *fGBTargets;  ///< Array of AliSurveyPoint of Button Targets
  TObjArray *fLBTargets;  ///< Array of TVector3 or AliSurveyPoint of local position of Button Targets
  TGeoCombiTrans *fLocalTrf; ///< Local transformation
  TGeoCombiTrans *fAlignTrf; ///< Local alignment transformation
  TGeoCombiTrans *fBaseTrf;  ///< Base Transformation

  Bool_t fOwnerLocalTrf;    ///< Flag for owner of fLocalTrf
  Bool_t fOwnerAlignTrf;    ///< Flag for owner of fAlignTrf
  Bool_t fOwnerBaseTrf;     ///< Flag for owner of fBaseTrf

  Bool_t fUseCM;            ///< Use centimeters, survey units are mm but aliroot uses cm

  TF2 *fPlane;  ///< TF2 for plane fitting

  TFitter *fFitter;  ///< Fitter for best local to global transformation

  Double_t fXMin;   
  Double_t fXMax;   
  Double_t fYMin;   
  Double_t fYMax;   
  Double_t fZMin;   
  Double_t fZMax;   


ClassDef(AliMUONSurveyObj, 0) //Class for alignment of muon spectrometer
};

#endif
