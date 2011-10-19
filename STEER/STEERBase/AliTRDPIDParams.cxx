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

#include "AliTRDPIDParams.h"

ClassImp(AliTRDPIDParams)

const Double_t AliTRDPIDParams::kVerySmall = 1e-5;

//____________________________________________________________
AliTRDPIDParams::AliTRDPIDParams():
  TNamed(),
  fEntries(NULL)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliTRDPIDParams::AliTRDPIDParams(const char *name) :
  TNamed(name, ""),
  fEntries(NULL)
{
  //
  // Default constructor
  //
  fEntries = new TSortedList;
}

//____________________________________________________________
AliTRDPIDParams::~AliTRDPIDParams(){
  //
  // Destructor
  //
  delete fEntries;
}

//____________________________________________________________
Bool_t AliTRDPIDParams::GetThresholdParameters(Int_t ntracklets, Double_t efficiency, Double_t *params) const{
  //
  // Retrieve params
  // Use IsEqual definition
  //
  AliTRDPIDThresholds test(ntracklets, efficiency);
  TObject *result = fEntries->FindObject(&test);
  if(!result){ 
    AliDebug(1, Form("No threshold params found for %d tracklets and an electron efficiency of %f", ntracklets, efficiency));
    return kFALSE;
  }
  AliTRDPIDThresholds *parResult = static_cast<AliTRDPIDThresholds *>(result);
  AliDebug(1, Form("Threshold params found: NTracklets %d, Electron Efficiency %f", parResult->GetNTracklets(), parResult->GetElectronEfficiency()));
  memcpy(params, parResult->GetThresholdParams(), sizeof(Double_t) * 4);
  return kTRUE;
}

//____________________________________________________________
void AliTRDPIDParams::SetThresholdParameters(Int_t ntracklets, Double_t effMin, Double_t effMax, Double_t *params){
  // 
  // Store new Params in the Object
  //
  if(effMin > effMax){
    AliError("Min. efficiency has to be >= max. efficiency");
    return;
  }
  AliDebug(1, Form("Save Parameters for %d tracklets at and electron efficiency of [%f|%f]", ntracklets, effMin, effMax));
  fEntries->Add(new AliTRDPIDThresholds(ntracklets, effMin, effMax, params));
}

//____________________________________________________________
void AliTRDPIDParams::Print(Option_t *) const {
  printf("Available thresholds:\n");
  printf("_________________________________________\n");
  TIter objects(fEntries);
  AliTRDPIDThresholds *par;
  while((par = dynamic_cast<AliTRDPIDThresholds *>(objects()))){
    printf("Number of tracklets %d, Electron efficiency %f\n", par->GetNTracklets(), par->GetElectronEfficiency());
  }
}

//____________________________________________________________
AliTRDPIDParams::AliTRDPIDThresholds::AliTRDPIDThresholds():
  TObject(),
  fNTracklets(0)
{
   // 
   // Default constructor
   //
   memset(fParams, 0, sizeof(Double_t) * 4);
   memset(fEfficiency, 0, sizeof(Double_t) * 2);
}

//____________________________________________________________
AliTRDPIDParams::AliTRDPIDThresholds::AliTRDPIDThresholds(Int_t nTracklets, Double_t effMin, Double_t effMax, Double_t *params) :
  TObject(),
  fNTracklets(nTracklets)
{
  //
  // Default Constructor 
  //
  fEfficiency[0] = effMin;
  fEfficiency[1] = effMax;  
  if(params) memcpy(fParams, params, sizeof(Double_t) * 4);
  else memset(fParams, 0, sizeof(Double_t) * 4);
}

//____________________________________________________________
AliTRDPIDParams::AliTRDPIDThresholds::AliTRDPIDThresholds(Int_t nTracklets, Double_t eff, Double_t *params) :
  TObject(),
  fNTracklets(nTracklets)
{
  //
  // Constructor used to find object in sorted list
  //
  fEfficiency[0] = fEfficiency[1] = eff;  
  if(params) memcpy(fParams, params, sizeof(Double_t) * 4);
  else memset(fParams, 0, sizeof(Double_t) * 4);
}

//____________________________________________________________
AliTRDPIDParams::AliTRDPIDThresholds::AliTRDPIDThresholds(const AliTRDPIDThresholds &ref) :
  TObject(ref),
  fNTracklets(ref.fNTracklets)
{
  //
  // Copy constructor
  //
  memcpy(fParams, ref.fParams, sizeof(Double_t) * 4);
  memcpy(fEfficiency, ref.fEfficiency, sizeof(Double_t) * 2);
}

//____________________________________________________________
AliTRDPIDParams::AliTRDPIDThresholds &AliTRDPIDParams::AliTRDPIDThresholds::operator=(const AliTRDPIDThresholds &ref){
  //
  // Assignment operator
  //
  TObject::operator=(ref);

  fNTracklets = ref.fNTracklets;
  memcpy(fEfficiency, ref.fEfficiency, sizeof(Double_t) * 2);
  memcpy(fParams, ref.fParams, sizeof(Double_t) * 4);
  return *this;
}
        
//____________________________________________________________
Int_t AliTRDPIDParams::AliTRDPIDThresholds::Compare(const TObject *ref) const{
  //
  // Compares two objects
  // Order:
  //   First compare number of tracklets, if they are equal compare electron efficiency
  //
  const AliTRDPIDThresholds *refObj = static_cast<const AliTRDPIDThresholds *>(ref);
  if(fNTracklets < refObj->GetNTracklets()) return -1;
  else if(fNTracklets > refObj->GetNTracklets()) return 1;
  else{
    if(fEfficiency[1] < refObj->GetElectronEfficiency(0)) return -1;
    else if(fEfficiency[0] > refObj->GetElectronEfficiency(1)) return 1;
    else return 0;
  }
}

//____________________________________________________________
Bool_t AliTRDPIDParams::AliTRDPIDThresholds::IsEqual(const TObject *ref) const {
  //
  // Check for equality 
  // Tracklets and Efficiency are used
  //
  const AliTRDPIDThresholds *refObj = dynamic_cast<const AliTRDPIDThresholds *>(ref);
  if(!refObj) return kFALSE;
  Bool_t eqNTracklets = fNTracklets == refObj->GetNTracklets();
  Bool_t eqEff = kFALSE;
  Bool_t hasRange = TMath::Abs(fEfficiency[1] - fEfficiency[0]) > kVerySmall;
  Bool_t hasRangeRef = TMath::Abs(refObj->GetElectronEfficiency(1) - refObj->GetElectronEfficiency(0)) > kVerySmall;
  if(hasRange && hasRangeRef){
    // Both have ranges, check if they match
    eqEff = TMath::Abs(fEfficiency[0] - refObj->GetElectronEfficiency(0)) < kVerySmall && TMath::Abs(fEfficiency[1] - refObj->GetElectronEfficiency(1)) < kVerySmall;
  } else if(hasRange){
    // this object has ranges, check if the efficiency of ref is inside the range
    eqEff = refObj->GetElectronEfficiency(0) >= fEfficiency[0] && refObj->GetElectronEfficiency(0) < fEfficiency[1];
  } else {
    // ref has ranges, check if this is in range
    eqEff = fEfficiency[0] >= refObj->GetElectronEfficiency(0) && fEfficiency[0] < refObj->GetElectronEfficiency(1);
  }
  
  return  eqNTracklets && eqEff;
}
