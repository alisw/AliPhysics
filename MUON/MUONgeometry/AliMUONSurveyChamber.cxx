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
/// \class AliMUONSurveyChamber
/// Class for the survey processing of the ALICE DiMuon spectrometer 
///
/// This object provides the methods specific to the chambers (frames)
///
/// \author Javier Castillo
//-----------------------------------------------------------------------------


#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TH2.h>

#include "AliLog.h"
#include "AliSurveyObj.h"
#include "AliSurveyPoint.h"

#include "AliMUONSurveyChamber.h"
#include "AliMUONSurveyDetElem.h"

/// \cond CLASSIMP
ClassImp(AliMUONSurveyChamber)
/// \endcond

AliMUONSurveyChamber::AliMUONSurveyChamber(Int_t lChamberId) 
  : AliMUONSurveyObj() /// Constructor with chamber id
  , fChamberId(lChamberId)
  , fNDetElem(0)
  , fSurveyObj(0x0)
  , fSurveyDetElem(0x0)
{
  /// Constructor with the chamber id
  fSurveyObj = new AliSurveyObj();
  fSurveyDetElem = new TClonesArray("AliMUONSurveyDetElem",4);
}

AliMUONSurveyChamber::~AliMUONSurveyChamber() {
  /// Destructor
  if (fSurveyObj) delete fSurveyObj;
  if (fSurveyDetElem) fSurveyDetElem->Delete();
}

Int_t AliMUONSurveyChamber::AddSurveyDetElem(Int_t lDetElemId) {
  /// Add a surveyed detection element to this chamber
  if (!fSurveyDetElem) {
    fSurveyDetElem = new TClonesArray("AliMUONSurveyDetElem",4);
    fNDetElem=0;
  }
  new((*fSurveyDetElem)[fNDetElem++]) AliMUONSurveyDetElem(lDetElemId,this);

  return fNDetElem;
}

AliMUONSurveyDetElem* AliMUONSurveyChamber::GetDetElem(Int_t lDetElemIndex) {
  /// Return AluMUONSurveyDetElem at lDetElemIndex
  if (lDetElemIndex<0||(lDetElemIndex>=fNDetElem)) {
    return 0x0;
  }
  else {
    return (AliMUONSurveyDetElem*)fSurveyDetElem->At(lDetElemIndex);
  }
}

Int_t AliMUONSurveyChamber::AddStickerTargets(TObjArray *pArray, TString stBaseName, Int_t lTargetMax) {
  /// Add a maximum of lTargetMax sticker targets with stBaseName from the pArray of targets
  return AliMUONSurveyObj::AddStickerTargets(pArray, stBaseName, lTargetMax);
}

Int_t AliMUONSurveyChamber::AddStickerTargets(TString stBaseName, Int_t lTargetMax) {
  /// Add a maximum of lTargetMax sticker targets with stBaseName from internal SurveyObj
  if (!fSurveyObj || !fSurveyObj->GetData()) {
    AliError("Survey data is missing!");
    return 0;    
  }
  return AddStickerTargets(fSurveyObj->GetData(),stBaseName,lTargetMax);
}

Int_t AliMUONSurveyChamber::AddGButtonTargets(TObjArray *pArray, TString stBaseName, Int_t lTargetMax) {
  /// Add a maximum of lTargetMax global targets with stBaseName from the pArray of targets
  return AliMUONSurveyObj::AddGButtonTargets(pArray, stBaseName, lTargetMax);
}

Int_t AliMUONSurveyChamber::AddGButtonTargets(TString stBaseName, Int_t lTargetMax) {
  /// Add a maximum of lTargetMax global button targets with stBaseName from internal SurveyObj
  if (!fSurveyObj || !fSurveyObj->GetData()) {
    AliError("Survey data is missing!");
    return 0;    
  }
  return AddGButtonTargets(fSurveyObj->GetData(),stBaseName,lTargetMax);
}

void AliMUONSurveyChamber::SetLocalTransformation(TGeoCombiTrans *localTrf, Bool_t ownerLocalTrf) {
  AliMUONSurveyObj::SetLocalTransformation(localTrf,ownerLocalTrf);
  AliMUONSurveyObj::SetBaseTransformation(localTrf,ownerLocalTrf);
}

void AliMUONSurveyChamber::PrintSurveyReport() {
  /// Print the survey report information and data
  printf("--> %d\n", fSurveyObj->GetEntries());

  printf("Title: \"%s\"\n", fSurveyObj->GetReportTitle().Data());
  printf("Date: \"%s\"\n", fSurveyObj->GetReportDate().Data());
  printf("Detector: \"%s\"\n", fSurveyObj->GetDetector().Data());
  printf("URL: \"%s\"\n", fSurveyObj->GetURL().Data());
  printf("Number: \"%d\"\n", fSurveyObj->GetReportNumber());
  printf("Version: \"%d\"\n", fSurveyObj->GetReportVersion());
  printf("Observations: \"%s\"\n", fSurveyObj->GetObservations().Data());
  printf("Coordinate System: \"%s\"\n", fSurveyObj->GetCoordSys().Data());
  printf("Measurement Units: \"%s\"\n", fSurveyObj->GetUnits().Data());
  printf("Nr Columns: \"%d\"\n", fSurveyObj->GetNrColumns());

  TObjArray *colNames = fSurveyObj->GetColumnNames();
  for (Int_t i = 0; i < colNames->GetEntries(); ++i)
    printf("  Column %d --> \"%s\"\n", i, ((TObjString *) colNames->At(i))->GetString().Data());

  // Get Array of surveyed points
  printf("Points:\n");
  TObjArray *points = fSurveyObj->GetData();
  
  for (Int_t i = 0; i < points->GetEntries(); ++i)
    printf("  Point %d --> \"%s\"  %s \n", i, ((AliSurveyPoint *) points->At(i))->GetPointName().Data(), points->At(i)->GetName());

}

void AliMUONSurveyChamber::FillCPSTHistograms(TString baseNameC, TH2 *hCPSTc, TString baseNameA, TH2 *hCPSTa) {
  /// Fill Chamber Plane Sticker Targest histograms for monitoring
  if(baseNameC.IsNull()||!hCPSTc){
    AliError("Need base name for points on side C and/or a histogram for them!");
    return;
  }
  AliMUONSurveyObj::FillSTHistograms(baseNameC,hCPSTc,baseNameA,hCPSTa);
}

void AliMUONSurveyChamber::FillDESTHistograms(TString baseNameC, TH2 *hDESTc, TString baseNameA, TH2 *hDESTa) {
  /// Fill Detection Element Sticker Targest histograms for monitoring
  for (Int_t iDE=0; iDE<GetNDetElem(); iDE++){
    GetDetElem(iDE)->FillSTHistograms(baseNameC,hDESTc,baseNameA,hDESTa);
  }
}

Double_t AliMUONSurveyChamber::GetMeanDetElemAlignResX() {
  /// Return the average uncertainty of the det. elem. translations along x parameter 
  Double_t alignResX = 0.;
  for (int iDE=0; iDE<GetNDetElem(); iDE++){
    alignResX += GetDetElem(iDE)->GetAlignResX();
  }
  if (GetNDetElem()==0){
    AliError("This Chamber has 0 detection elements!");
    return 0.;
  }
  return alignResX/GetNDetElem();
}

Double_t AliMUONSurveyChamber::GetMeanDetElemAlignResY() {
  /// Return the average uncertainty of the det. elem. translations along y parameter 
  Double_t alignResY = 0.;
  for (int iDE=0; iDE<GetNDetElem(); iDE++){
    alignResY += GetDetElem(iDE)->GetAlignResY();
  }
  if (GetNDetElem()==0){
    AliError("This Chamber has 0 detection elements!");
    return 0.;
  }
  return alignResY/GetNDetElem();
}
