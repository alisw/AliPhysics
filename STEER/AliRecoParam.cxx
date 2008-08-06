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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// ALICE Reconstruction parameterization:                                    //
//                                                                           //
//                                                                           //
// Base Class for Detector reconstruction parameters                         //
// Revision: cvetan.cheshkov@cern.ch 12/06/2008                              //
// Its structure has been revised and it is interfaced to AliEventInfo.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"
#include "AliDetectorRecoParam.h"

#include "AliLog.h"
#include "AliRecoParam.h"

ClassImp(AliRecoParam)

AliRecoParam::AliRecoParam(): 
  TObject(),
  fEventSpecie(kDefault)
{
  // Default constructor
  // ...
  for(Int_t iDet = 0; iDet < kNDetectors; iDet++)
    fDetRecoParams[iDet] = NULL;
  for(Int_t iSpecie = 0; iSpecie < kNSpecies; iSpecie++) {
    for(Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      fDetRecoParamsIndex[iSpecie][iDet] = -1;
    }
  }
}

AliRecoParam::AliRecoParam(const AliRecoParam& par) :
  TObject(),
  fEventSpecie(par.fEventSpecie)
{
  // copy constructor
  for(Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (par.fDetRecoParams[iDet])
      fDetRecoParams[iDet] = (TObjArray*)(par.fDetRecoParams[iDet]->Clone());
    else
      fDetRecoParams[iDet] = NULL;
  }
  for(Int_t iSpecie = 0; iSpecie < kNSpecies; iSpecie++) {
    for(Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      fDetRecoParamsIndex[iSpecie][iDet] = par.fDetRecoParamsIndex[iSpecie][iDet];
    }
  }
}

//_____________________________________________________________________________
AliRecoParam& AliRecoParam::operator = (const AliRecoParam& par)
{
  // assignment operator

  if(&par == this) return *this;

  this->~AliRecoParam();
  new(this) AliRecoParam(par);
  return *this;
}

AliRecoParam::~AliRecoParam(){
  // Destructor
  // ...
  // Delete the array with the reco-param objects
  for(Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (fDetRecoParams[iDet]){
      fDetRecoParams[iDet]->Delete();
      delete fDetRecoParams[iDet];
    }
  }
}

void  AliRecoParam::Print(Option_t *option) const {
  //
  // Print reconstruction setup
  //
  for(Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (fDetRecoParams[iDet]){
      printf("AliDetectorRecoParam objects for detector %d:\n",iDet); 
      Int_t nparam = fDetRecoParams[iDet]->GetEntriesFast();
      for (Int_t iparam=0; iparam<nparam; iparam++){
	AliDetectorRecoParam * param = (AliDetectorRecoParam *)fDetRecoParams[iDet]->At(iparam);
	if (!param) continue;
	param->Print(option);
      }
    }
    else {
      printf("No AliDetectorRecoParam objects specified for detector %d\n",iDet); 
    }
  }
}

void AliRecoParam::SetEventSpecie(const AliRunInfo */*runInfo*/, const AliEventInfo &/*evInfo*/)
{
  // To be implemented
  // Here we return always kDefault!!
  fEventSpecie = kDefault;
}

const AliDetectorRecoParam *AliRecoParam::GetDetRecoParam(Int_t iDet) const
{
  // Return AliDetectorRecoParam object for a given detector
  // according to the event specie provided as an argument
  if (!fDetRecoParams[iDet]) return NULL;
  if (fDetRecoParams[iDet]->GetEntries() == 0) return NULL;

  for(Int_t iBit = 0; iBit < kNSpecies; iBit++) {
    if (fEventSpecie & (1 << iBit)) {
      if (fDetRecoParamsIndex[iBit][iDet] >= 0)
	return (AliDetectorRecoParam *)fDetRecoParams[iDet]->At(fDetRecoParamsIndex[iBit][iDet]);
      else
	return (AliDetectorRecoParam *)fDetRecoParams[iDet]->At(fDetRecoParamsIndex[0][iDet]);
    }
  }

  // Default one
  AliError(Form("Invalid event specie: %d!",fEventSpecie));
  return (AliDetectorRecoParam *)fDetRecoParams[iDet]->At(fDetRecoParamsIndex[0][iDet]);
}

void  AliRecoParam::AddDetRecoParam(Int_t iDet, AliDetectorRecoParam* param)
{
  // Add an instance of reco params object into
  // the fDetRecoParams for detector iDet
  // Updates the fDetRecoParams index
  if (!fDetRecoParams[iDet]) fDetRecoParams[iDet] = new TObjArray;
  fDetRecoParams[iDet]->AddLast(param);
  Int_t index = fDetRecoParams[iDet]->GetLast();

  // Index
  Int_t specie = param->GetEventSpecie();
  for(Int_t iBit = 0; iBit < kNSpecies; iBit++) {
    if (specie & (1 << iBit)) {
      fDetRecoParamsIndex[iBit][iDet] = index;
    }
  }
}

Bool_t AliRecoParam::AddDetRecoParamArray(Int_t iDet, TObjArray* parArray)
{
  // Add an array of reconstruction parameter objects
  // for a given detector
  // Basic check on the consistency of the array
  Bool_t defaultFound = kFALSE;
  for(Int_t i = 0; i < parArray->GetEntriesFast(); i++) {
    AliDetectorRecoParam *par = (AliDetectorRecoParam*)parArray->At(i);
    if (!par) continue;
    if (par->IsDefault()) defaultFound = kTRUE;

    Int_t specie = par->GetEventSpecie();
    for(Int_t iBit = 0; iBit < kNSpecies; iBit++) {
      if (specie & (1 << iBit)) {
	fDetRecoParamsIndex[iBit][iDet] = i;
      }
    }
 }
   
  fDetRecoParams[iDet] = parArray;

  return defaultFound;
}
