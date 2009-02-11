/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-----------------------------------------------------------------------------
/// \class AliMUONSurveyDetElem
/// Class for the survey processing of the ALICE DiMuon spectrometer 
///
/// This object provides the methods specific to the detection elements
///
/// \author Javier Castillo
//-----------------------------------------------------------------------------

#include "TGeoMatrix.h"

#include "AliLog.h"
#include "AliSurveyObj.h"

#include "AliMUONSurveyChamber.h"
#include "AliMUONSurveyDetElem.h"

/// \cond CLASSIMP
ClassImp(AliMUONSurveyDetElem)
/// \endcond

AliMUONSurveyDetElem::AliMUONSurveyDetElem(Int_t lDetElemId) 
  : AliMUONSurveyObj() 
  , fDetElemId(lDetElemId)
  , fSurveyChamber(0x0)
{
/// Constructor with detection element id
}

AliMUONSurveyDetElem::AliMUONSurveyDetElem(Int_t lDetElemId, AliMUONSurveyChamber *lSurveyChamber) 
  : AliMUONSurveyObj() /// Constructor with mother chamber provided
  , fDetElemId(lDetElemId)
  , fSurveyChamber(lSurveyChamber)
{
/// Constructor with mother chamber provided
}

AliMUONSurveyDetElem::~AliMUONSurveyDetElem() {
  /// Destructor
}

Int_t AliMUONSurveyDetElem::AddStickerTargets(TObjArray *pArray, TString stBaseName, Int_t lTargetMax) {
  return AliMUONSurveyObj::AddStickerTargets(pArray, stBaseName, lTargetMax);
}

Int_t AliMUONSurveyDetElem::AddStickerTargets(TString stBaseName, Int_t lTargetMax) {
  /// Add a maximum of lTargetMax sticker targets with stBaseName from internal SurveyObj
  if (!fSurveyChamber) {
    AliError("Pointer to mother chamber has not been set!");
    return 0;
  }
  if (!fSurveyChamber->GetSurveyObj() || !fSurveyChamber->GetSurveyObj()->GetData()) {
    AliError("Survey data is missing!");
    return 0;    
  }
  return AddStickerTargets(fSurveyChamber->GetSurveyObj()->GetData(),stBaseName,lTargetMax);
}

Int_t AliMUONSurveyDetElem::AddGButtonTargets(TObjArray *pArray, TString stBaseName, Int_t lTargetMax) {
  return AliMUONSurveyObj::AddGButtonTargets(pArray, stBaseName, lTargetMax);
}

Int_t AliMUONSurveyDetElem::AddGButtonTargets(TString stBaseName, Int_t lTargetMax) {
  /// Add a maximum of lTargetMax sticker targets with stBaseName from internal SurveyObj
  if (!fSurveyChamber) {
    AliError("Pointer to mother chamber has not been set!");
    return 0;
  }
  if (!fSurveyChamber->GetSurveyObj() || !fSurveyChamber->GetSurveyObj()->GetData()) {
    AliError("Survey data is missing!");
    return 0;    
  }
  return AddGButtonTargets(fSurveyChamber->GetSurveyObj()->GetData(),stBaseName,lTargetMax);
}

void AliMUONSurveyDetElem::SetLocalTransformation(TGeoCombiTrans *localTrf, Bool_t ownerLocalTrf) {
  /// Set the geometry transformation of this detection element
  AliMUONSurveyObj::SetLocalTransformation(localTrf,ownerLocalTrf);
  if (!fSurveyChamber) {
    AliWarning("Pointer to mother chamber has not been set!");
    AliMUONSurveyObj::SetBaseTransformation(localTrf,ownerLocalTrf);
  } else {
    if (fSurveyChamber->GetLocalTrf()){
      if (fSurveyChamber->GetAlignTrf()){
	AliMUONSurveyObj::SetBaseTransformation(new TGeoCombiTrans((*(fSurveyChamber->GetLocalTrf()))*(*(fSurveyChamber->GetAlignTrf()))*(*localTrf)),kTRUE);
      } else {
	AliWarning("Mother chamber has not been aligned yet!");
	AliMUONSurveyObj::SetBaseTransformation(new TGeoCombiTrans((*(fSurveyChamber->GetLocalTrf()))*(*localTrf)),kTRUE);
      }
    } else {
      AliWarning("Mother chamber has no local transformation");
      AliMUONSurveyObj::SetBaseTransformation(localTrf,ownerLocalTrf);      
    }
  }
}

void AliMUONSurveyDetElem::PrintLocalTrf() {
  printf("DetElem%d Th",fDetElemId);
  AliMUONSurveyObj::PrintLocalTrf();
}

void AliMUONSurveyDetElem::PrintAlignTrf() {
  printf("DetElem%d d",fDetElemId);
  AliMUONSurveyObj::PrintAlignTrf();
}
