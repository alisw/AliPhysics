#ifndef ALIMUONSURVEYDETELEM_H
#define ALIMUONSURVEYDETELEM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup geometry
/// \class AliMUONSurveyDetElem
/// \brief Class for survey of detection elements of the muon spectrometer
//
// Author: Javier Castillo

#include "AliMUONSurveyObj.h"

class AliMUONSurveyChamber;

class AliMUONSurveyDetElem:public AliMUONSurveyObj
{

 public:
  AliMUONSurveyDetElem(Int_t lDetElemId);
  AliMUONSurveyDetElem(Int_t lDetElemId, AliMUONSurveyChamber *lSurveyChamber);

  virtual Int_t AddStickerTargets(TString stBaseName, Int_t lTargetMax = 9);
  virtual Int_t AddGButtonTargets(TString btBaseName, Int_t lTargetMax = 9);

  virtual Int_t AddStickerTargets(TObjArray *pArray, TString stBaseName, Int_t lTargetMax = 9);
  virtual Int_t AddGButtonTargets(TObjArray *pArray, TString btBaseName, Int_t lTargetMax = 9);

  virtual ~AliMUONSurveyDetElem();
 
  virtual void SetLocalTransformation(TGeoCombiTrans *localTrf, Bool_t ownerLocalTrf = kFALSE);

  virtual void PrintLocalTrf();
  virtual void PrintAlignTrf();

 private:
  /// Not implemented
  AliMUONSurveyDetElem(const AliMUONSurveyDetElem& right);
  /// Not implemented
  AliMUONSurveyDetElem&  operator = (const AliMUONSurveyDetElem& right);

  Int_t fDetElemId;   ///< Detection element id
  AliMUONSurveyChamber *fSurveyChamber;  ///< Pointer to mother survey chamber object


ClassDef(AliMUONSurveyDetElem, 0) //Class for survey det. elem. of muon spectrometer
};

#endif
