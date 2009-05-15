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

//-------------------------------------------------------------------------
//                      Class AliRsnPIDDefESD
//                     -------------------
//           Simple collection of reconstructed tracks
//           selected from an ESD event
//           to be used for analysis.
//           .........................................
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#include <TMath.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliRsnDaughter.h"
#include "AliRsnPIDDefESD.h"

ClassImp(AliRsnPIDDefESD)

//_____________________________________________________________________________
AliRsnPIDDefESD::AliRsnPIDDefESD() : fUseESDWeights(kTRUE) {
//
// Default constructor.
// By default, it is set for using ESD weights,
// so, values in other members are meaningless.
//

  Int_t i;
  for (i = 0; i < kDetectors; i++) {
    fUseDet[i] = kTRUE;
    fDivValue[i] = 0.0;
    fUseHigher[i] = kTRUE;
  }
}

//_____________________________________________________________________________
AliRsnPIDDefESD::AliRsnPIDDefESD(const AliRsnPIDDefESD& copy) :
    TObject(copy),
    fUseESDWeights(copy.fUseESDWeights) {
//
// Copy constructor.
// Implemented to manage passing of this object to functions
//

  Int_t i;
  for (i = 0; i < kDetectors; i++) {
    fUseDet[i] = copy.fUseDet[i];
    fDivValue[i] = copy.fDivValue[i];
    fUseHigher[i] = copy.fUseHigher[i];
  }
}

//_____________________________________________________________________________
void AliRsnPIDDefESD::SetScheme(EScheme scheme, Double_t divValue) {
//
// Set one of the predefined schemes
//

  switch (scheme) {
    case kSchemeESD:
      fUseESDWeights = kTRUE;
      break;
    case kSchemeITS:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kITS);
      SetDivValue(kITS, 0.0);
      break;
    case kSchemeTPC:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kTPC);
      SetDivValue(kTPC, 0.0);
      break;
    case kSchemeTOF:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kTOF);
      SetDivValue(kTOF, 0.0);
      break;
    case kSchemeITSandTPC:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kITS);
      IncludeDet(kTPC);
      SetDivValue(kITS, 0.0);
      SetDivValue(kTPC, 0.0);
      break;
    case kSchemeITSandTOF:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kITS);
      IncludeDet(kTOF);
      SetDivValue(kITS, 0.0);
      SetDivValue(kTOF, 0.0);
      break;
    case kSchemeTPCandTOF:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kTPC);
      IncludeDet(kTOF);
      SetDivValue(kTPC, 0.0);
      SetDivValue(kTOF, 0.0);
      break;
    case kSchemeITSandTPCandTOF:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kITS);
      IncludeDet(kTPC);
      IncludeDet(kTOF);
      SetDivValue(kITS, 0.0);
      SetDivValue(kTPC, 0.0);
      SetDivValue(kTOF, 0.0);
      break;
    case kSchemeITSandTPCandTOFwithSP:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kITS);
      IncludeDet(kTPC);
      IncludeDet(kTOF);
      SetDivValue(kITS, 0.0);
      SetDivValue(kTPC, 0.0);
      SetDivValue(kTOF, divValue);
      break;
    case kSchemeITSandTPCorTOFwithSP:
      fUseESDWeights = kFALSE;
      ExcludeAll();
      IncludeDet(kITS);
      IncludeDet(kTPC);
      IncludeDet(kTOF);
      SetDivValue(kITS, divValue, kFALSE);
      SetDivValue(kTPC, divValue, kFALSE);
      SetDivValue(kTOF, divValue, kTRUE);
      break;
    default:
      AliWarning("PID scheme unrecognized. Set to ESD");
      fUseESDWeights = kTRUE;
  }
}

//_____________________________________________________________________________
void AliRsnPIDDefESD::ComputeWeights(AliESDtrack *track, Double_t *weights) {
//
// Computes the global PID weights using the given ranges
//

  if (fUseESDWeights) {
    track->GetESDpid(weights);
    return;
  }

  Double_t pt = track->Pt();
  Double_t w[kDetectors][AliPID::kSPECIES];
  track->GetITSpid(w[kITS]);
  track->GetTPCpid(w[kTPC]);
  track->GetTRDpid(w[kTRD]);
  track->GetTOFpid(w[kTOF]);
  track->GetHMPIDpid(w[kHMPID]);

  Int_t i, j;
  for (i = 0; i < kDetectors; i++) {
//     if (!fUseDet[i] || pt < fDivValue[i])
    if (!fUseDet[i] || !CheckDivValue((EDetector)i,pt)) {
      for (j = 0; j < AliPID::kSPECIES; j++) {
        w[i][j] = 1.0;
      }
    }
  }

  for (i = 0; i < AliPID::kSPECIES; i++) {
    weights[i] = w[kITS][i] * w[kTPC][i] * w[kTRD][i] * w[kTOF][i] * w[kHMPID][i];
  }
}

//_____________________________________________________________________________
void AliRsnPIDDefESD::PrintStatus() {
//
// Print informations about this object configurations
//

  AliInfo("===== PIDDef status messages -- BEGIN");

  if (fUseESDWeights) {
    AliInfo("Using ESD weights");
  } else {
    AliInfo("NOT using ESD weights");
  }

  Int_t i;
  for (i = 0; i < kDetectors; i++) {
    AliInfo(Form("Detector name: %s -- accepted: %s -- divValue = %3.1f useHigher = %s", DetName((EDetector)i), (fUseDet[i]?"YES":"NO"), fDivValue[i],(fUseHigher[i]?"YES":"NO")));
  }

  AliInfo("===== PIDDef status messages -- END");
}

//_____________________________________________________________________________
const char* AliRsnPIDDefESD::DetName(EDetector det) {
//
// Detector name for messages
//

  switch (det) {
    case kITS: return "ITS";
    case kTPC: return "TPC";
    case kTRD: return "TRD";
    case kTOF: return "TOF";
    case kHMPID: return "HMPID";
    default: return "undef";
  }
}

//_____________________________________________________________________________
void AliRsnPIDDefESD::SetDivValue(EDetector det, Double_t value, Bool_t userHigher) {
//
// Sets div.value properties for detector
//
  if (CheckBounds(det)) {
    fDivValue[det] = value;
    fUseHigher[det] = userHigher;
  }
}

//_____________________________________________________________________________
Bool_t AliRsnPIDDefESD::CheckDivValue(EDetector det,Double_t value) {
//
// Sets div.value properties for detector
//
  if (CheckBounds(det)) {
    if (fUseHigher[det]) {
      if (value > fDivValue[det]) return kTRUE;
      else return kFALSE;
    } else {
      if (value < fDivValue[det]) return kTRUE;
      else return kFALSE;
    }
  }

  return kTRUE;
}
