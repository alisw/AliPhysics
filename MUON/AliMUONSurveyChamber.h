#ifndef ALIMUONSURVEYCHAMBER_H
#define ALIMUONSURVEYCHAMBER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup geometry
/// \class AliMUONSurveyChamber
/// \brief Class for survey of chambers (frames) of the muon spectrometer
//
// Authors: Javier Castillo

#include "AliMUONSurveyObj.h"

class TClonesArray;
class TH2;

class AliSurveyObj;

class AliMUONSurveyDetElem;

class AliMUONSurveyChamber: public AliMUONSurveyObj
{

 public:
  AliMUONSurveyChamber(Int_t lChamberId);
  virtual ~AliMUONSurveyChamber();

  virtual Int_t AddStickerTargets(TString stBaseName, Int_t lTargetMax = 9);
  virtual Int_t AddGButtonTargets(TString btBaseName, Int_t lTargetMax = 9);

  virtual Int_t AddStickerTargets(TObjArray *pArray, TString stBaseName, Int_t lTargetMax = 9);
  virtual Int_t AddGButtonTargets(TObjArray *pArray, TString btBaseName, Int_t lTargetMax = 9);
 
  Int_t AddSurveyDetElem(Int_t lDetElemId);
  Int_t GetNDetElem() const  {return fNDetElem;}
  AliMUONSurveyDetElem* GetDetElem(Int_t lDetElemIndex);

  AliSurveyObj* GetSurveyObj() const {return fSurveyObj;}

  virtual void SetLocalTransformation(TGeoCombiTrans *localTrf, Bool_t ownerLocalTrf = kFALSE);

  void PrintSurveyReport();

  void FillCPSTHistograms(TString baseNameC, TH2 *hCPSTc, TString baseNameA="", TH2 *hCPSTa = 0);
  void FillDESTHistograms(TString baseNameC, TH2 *hCPSTc, TString baseNameA="", TH2 *hCPSTa = 0);

  Double_t GetMeanDetElemAlignResX(); 
  Double_t GetMeanDetElemAlignResY(); 

 private:
  /// Not implemented
  AliMUONSurveyChamber(const AliMUONSurveyChamber& right);
  /// Not implemented
  AliMUONSurveyChamber&  operator = (const AliMUONSurveyChamber& right);

  Int_t fChamberId;    ///< Chamber Id
  Int_t fNDetElem;     ///< Number of detection elements

  AliSurveyObj *fSurveyObj;     ///< Survey object containing the measurment
  TClonesArray *fSurveyDetElem; ///< Array of AliMUONSurveyDetElem

ClassDef(AliMUONSurveyChamber, 0) //Class for survey of muon spectrometer chambers
};

#endif
