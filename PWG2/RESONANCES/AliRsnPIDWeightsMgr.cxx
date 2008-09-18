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
//                      Class AliRsnPIDWeightsMgr
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
#include "AliRsnPIDWeightsMgr.h"

ClassImp(AliRsnPIDWeightsMgr)

//_____________________________________________________________________________
AliRsnPIDWeightsMgr::AliRsnPIDWeightsMgr()
{
//
// Default and unique constructor which initializes arrays
//
  Int_t i, j;
  for (i = 0; i < kDetectors; i++)
  {
    for (j = 0; j < AliRsnPID::kSpecies; j++) fWeights[i][j] = 0.0;
    fDetPtMin[i] = 0.0;
    fDetPtMax[i] = 1E25;
    fUseDet[i] = kTRUE;
  }
}

//_____________________________________________________________________________
void AliRsnPIDWeightsMgr::SetDetectorWeights(EDetector det, Double_t *weights)
{
//
// Copy into the array the PID weights of a detector
//

  Int_t i;
  if (!CheckBounds(det)) return;
  for (i = 0; i < AliRsnPID::kSpecies; i++) fWeights[det][i] = weights[i];
}

//_____________________________________________________________________________
void AliRsnPIDWeightsMgr::SetAcceptanceRange(EDetector det, Double_t ptmin, Double_t ptmax)
{
//
// Sets a range in pt for which the PID weights of this detector
// are accepted in the global computation, in the case that one
// did not accept to use the global ESD pid.
//

  if (!CheckBounds(det)) return;

  // swap values if necessary
  if (ptmin > ptmax)
  {
    AliWarning("Values passed in wrong order. Swapping");
    Double_t temp = ptmin;
    ptmin = ptmax;
    ptmax = temp;
  }
  fDetPtMin[det] = ptmin;
  fDetPtMax[det] = ptmax;
}

//_____________________________________________________________________________
Double_t AliRsnPIDWeightsMgr::GetWeight(AliRsnPID::EType type, Double_t pt)
{
//
// Computes the global PID weights using the given ranges
//

  if (type < 0 or type >= AliRsnPID::kSpecies)
  {
    AliError("Index out of range");
    return kFALSE;
  }

  Int_t i;
  Double_t prob = 1.0;
  for (i = 0; i < kDetectors; i++)
  {
    //AliInfo(Form("weights[%d] = %f %f %f %f %f", i, fWeights[i][0], fWeights[i][1], fWeights[i][2], fWeights[i][3], fWeights[i][4], fWeights[i][5]));
    if (!fUseDet[i]) continue;
    if (pt < fDetPtMin[i] || pt > fDetPtMax[i]) continue;
    prob *= fWeights[i][type];
  }

  return prob;
}

//_____________________________________________________________________________
Bool_t AliRsnPIDWeightsMgr::CheckBounds(EDetector det)
{
//
// Bounds checker
//

  if (det < 0 or det >= kDetectors)
  {
    AliError("Index out of range");
    return kFALSE;
  }

  return kTRUE;
}
