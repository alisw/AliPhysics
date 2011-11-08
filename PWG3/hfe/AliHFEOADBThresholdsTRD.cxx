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
//
// Container for TRD Threshold parameters stored in the OADB
//
// Author: Markus Fasel <M.Fasel@gsi.de>
//
#include <TMath.h>
#include <TSortedList.h>

#include "AliLog.h"

#include "AliHFEOADBThresholdsTRD.h"

ClassImp(AliHFEOADBThresholdsTRD)

const Double_t AliHFEOADBThresholdsTRD::fgkVerySmall = 1e-5;

//____________________________________________________________
AliHFEOADBThresholdsTRD::AliHFEOADBThresholdsTRD():
  TNamed(),
  fEntries(NULL)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliHFEOADBThresholdsTRD::AliHFEOADBThresholdsTRD(const char *name) :
  TNamed(name, ""),
  fEntries(NULL)
{
  //
  // Default constructor
  //
  fEntries = new TSortedList;
}

//____________________________________________________________
AliHFEOADBThresholdsTRD::~AliHFEOADBThresholdsTRD(){
  //
  // Destructor
  //
  delete fEntries;
}

//____________________________________________________________
Bool_t AliHFEOADBThresholdsTRD::GetThresholdParameters(Int_t ntracklets, Double_t efficiency, Double_t *params){
  //
  // Retrieve params
  // Use IsEqual definition
  //
  AliHFEthresholdParamsTRD test(ntracklets, efficiency);
  TObject *result = fEntries->FindObject(&test);
  if(!result){ 
    AliDebug(1, Form("No threshold params found for %d tracklets and an electron efficiency of %f", ntracklets, efficiency));
    return kFALSE;
  }
  AliHFEthresholdParamsTRD *parResult = static_cast<AliHFEthresholdParamsTRD *>(result);
  AliDebug(1, Form("Threshold params found: NTracklets %d, Electron Efficiency %f", parResult->GetNTracklets(), parResult->GetElectronEfficiency()));
  memcpy(params, parResult->GetThresholdParams(), sizeof(Double_t) * 4);
  return kTRUE;
}

//____________________________________________________________
void AliHFEOADBThresholdsTRD::SetThresholdParameters(Int_t ntracklets, Double_t efficiency, Double_t *params){
  // 
  // Store new Params in the Object
  //
  AliDebug(1, Form("Save Parameters for %d tracklets at and electron efficiency of %f", ntracklets, efficiency));
  fEntries->Add(new AliHFEthresholdParamsTRD(ntracklets, efficiency, params));
}

//____________________________________________________________
void AliHFEOADBThresholdsTRD::Print(Option_t *) const {
  //
  // Print thresholds
  //
  printf("Available thresholds:\n");
  printf("_________________________________________\n");
  TIter objects(fEntries);
  AliHFEthresholdParamsTRD *par;
  while((par = dynamic_cast<AliHFEthresholdParamsTRD *>(objects()))){
    printf("Number of tracklets %d, Electron efficiency %f\n", par->GetNTracklets(), par->GetElectronEfficiency());
  }
}

//____________________________________________________________
AliHFEOADBThresholdsTRD::AliHFEthresholdParamsTRD::AliHFEthresholdParamsTRD():
  TObject(),
  fNTracklets(0),
  fEfficiency(0.)
{
   // 
   // Default constructor
   //
   memset(fParams, 0, sizeof(Double_t) * 4);
}

//____________________________________________________________
AliHFEOADBThresholdsTRD::AliHFEthresholdParamsTRD::AliHFEthresholdParamsTRD(Int_t nTracklets, Double_t eff, Double_t *params) :
  TObject(),
  fNTracklets(nTracklets),
  fEfficiency(eff)
{
  //
  // Default Constructor
  //
  if(params) memcpy(fParams, params, sizeof(Double_t) * 4);
  else memset(fParams, 0, sizeof(Double_t) * 4);
}

//____________________________________________________________
AliHFEOADBThresholdsTRD::AliHFEthresholdParamsTRD::AliHFEthresholdParamsTRD(const AliHFEthresholdParamsTRD &ref) :
  TObject(ref),
  fNTracklets(ref.fNTracklets),
  fEfficiency(ref.fEfficiency)
{
  //
  // Copy constructor
  //
  memcpy(fParams, ref.fParams, sizeof(Double_t) * 4);
}

//____________________________________________________________
AliHFEOADBThresholdsTRD::AliHFEthresholdParamsTRD &AliHFEOADBThresholdsTRD::AliHFEthresholdParamsTRD::operator=(const AliHFEthresholdParamsTRD &ref){
  //
  // Assignment operator
  //
  TObject::operator=(ref);

  fNTracklets = ref.fNTracklets;
  fEfficiency = ref.fEfficiency;
  memcpy(fParams, ref.fParams, sizeof(Double_t) * 4);
  return *this;
}
        
//____________________________________________________________
Int_t AliHFEOADBThresholdsTRD::AliHFEthresholdParamsTRD::Compare(const TObject *ref) const{
  //
  // Compares two objects
  // Order:
  //   First compare number of tracklets, if they are equal compare electron efficiency
  //
  const AliHFEthresholdParamsTRD *refObj = static_cast<const AliHFEthresholdParamsTRD *>(ref);
  if(fNTracklets < refObj->GetNTracklets()) return -1;
  else if(fNTracklets > refObj->GetNTracklets()) return 1;
  else{
    if(fEfficiency < refObj->GetElectronEfficiency()) return -1;
    else if(fEfficiency > refObj->GetElectronEfficiency()) return 1;
    else return 0;
  }
}

//____________________________________________________________
Bool_t AliHFEOADBThresholdsTRD::AliHFEthresholdParamsTRD::IsEqual(const TObject *ref) const {
  //
  // Check for equality 
  // Tracklets and Efficiency are used
  //
  const AliHFEthresholdParamsTRD *refObj = dynamic_cast<const AliHFEthresholdParamsTRD *>(ref);
  if(!refObj) return kFALSE;
  return fNTracklets == refObj->GetNTracklets() && TMath::Abs(fEfficiency -  refObj->GetElectronEfficiency()) < fgkVerySmall;
}
