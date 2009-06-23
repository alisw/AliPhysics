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
/************************************************************************
 *                                                                      *
 * Class for TRD PID                                                    *
 * Implements the abstract base class AliHFEpidbase        *
 * Make PID does the PID decision                                       *
 * Class further contains TRD specific cuts and QA histograms           *
 *                                                                      *
 * Authors:                                                             *
 *   Markus Fasel <M.Fasel@gsi.de>                                      *
 *                                                                      *
 ************************************************************************/
#include <TAxis.h>
#include <TFile.h>
#include <TH1F.h>
#include <TIterator.h>
#include <TKey.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>

#include "AliESDtrack.h"
#include "AliPID.h"

#include "AliHFEpidTRD.h"

ClassImp(AliHFEpidTRD)

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const char* name) :
    AliHFEpidBase(name)
  , fPIDMethod(kNN)
{
  //
  // default  constructor
  // 
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
}

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const AliHFEpidTRD &ref):
    AliHFEpidBase("")
  , fPIDMethod(kLQ)
{
  //
  // Copy constructor
  //
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
  ref.Copy(*this);
}

//___________________________________________________________________
AliHFEpidTRD &AliHFEpidTRD::operator=(const AliHFEpidTRD &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    ref.Copy(*this);
  }
  return *this;
}

//___________________________________________________________________
void AliHFEpidTRD::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidTRD &target = dynamic_cast<AliHFEpidTRD &>(ref);

  target.fPIDMethod = fPIDMethod;
  memcpy(target.fThreshParams, fThreshParams, sizeof(Double_t) * kThreshParams);
  AliHFEpidBase::Copy(ref);
}

//___________________________________________________________________
AliHFEpidTRD::~AliHFEpidTRD(){
  //
  // Destructor
  //
}

//______________________________________________________
Bool_t AliHFEpidTRD::InitializePID(){
  //
  // InitializePID: Load TRD thresholds and create the electron efficiency axis
  // to navigate 
  //
  InitParameters();
  return kTRUE;
}

//______________________________________________________
Int_t AliHFEpidTRD::IsSelected(AliVParticle *track){
  //
  // Does PID for TRD alone:
  // PID thresholds based on 90% Electron Efficiency level approximated by a linear 
  // step function
  //
  AliESDtrack *esdTrack = 0x0;
  if(!(esdTrack = dynamic_cast<AliESDtrack *>(track))) return kFALSE;
  Double_t p = esdTrack->GetOuterParam() ? esdTrack->GetOuterParam()->P() : esdTrack->P();
  if(p < 0.6) return 0;

  Double_t pidProbs[AliPID::kSPECIES];
  esdTrack->GetTRDpid(pidProbs);
  if(pidProbs[AliPID::kElectron] > GetTRDthresholds(0.91, p)) return 11;
  return 0;
}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetTRDthresholds(Double_t electronEff, Double_t p){ 
  //
  // Return momentum dependent and electron efficiency dependent TRD thresholds
  // 
  Double_t params[4];
  GetParameters(electronEff, params);
  return 1. - params[0] * p - params[1] * TMath::Exp(-params[2] * p) - params[3];
}

//___________________________________________________________________
void AliHFEpidTRD::InitParameters(){
  //
  // Fill the Parameters into an array
  //

  // Parameters for 6 Layers
  fThreshParams[0] = -0.001839; // 0.7 electron eff
  fThreshParams[1] = 0.000276;
  fThreshParams[2] = 0.044902; 
  fThreshParams[3] = 1.726751;
  fThreshParams[4] = -0.002405; // 0.75 electron eff
  fThreshParams[5] = 0.000372;
  fThreshParams[6] = 0.061775;
  fThreshParams[7] = 1.739371;
  fThreshParams[8] = -0.003178; // 0.8 electron eff
  fThreshParams[9] = 0.000521;
  fThreshParams[10] = 0.087585;
  fThreshParams[11] = 1.749154;
  fThreshParams[12] = -0.004058; // 0.85 electron eff
  fThreshParams[13] = 0.000748;
  fThreshParams[14] = 0.129583;
  fThreshParams[15] = 1.782323;
  fThreshParams[16] = -0.004967; // 0.9 electron eff
  fThreshParams[17] = 0.001216;
  fThreshParams[18] = 0.210128;
  fThreshParams[19] = 1.807665;
  fThreshParams[20] = -0.000996; // 0.95 electron eff
  fThreshParams[21] = 0.002627;
  fThreshParams[22] = 0.409099;
  fThreshParams[23] = 1.787076;
}

//___________________________________________________________________
void AliHFEpidTRD::GetParameters(Double_t electronEff, Double_t *parameters){
  //
  // return parameter set for the given efficiency bin
  //
  Int_t effbin = static_cast<Int_t>((electronEff - 0.7)/0.3 * 6.);
  memcpy(parameters, fThreshParams + effbin * 4, sizeof(Double_t) * 4);
}
