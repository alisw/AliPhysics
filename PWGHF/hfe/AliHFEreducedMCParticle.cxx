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
// Debug MC particle
// the tree is represented as reduced events
// 
// Authors:
//   M.Fasel <M.Fasel@gsi.de>
//
#include "TObjArray.h"
#include <cstring>

#include "AliHFEreducedMCParticle.h"

ClassImp(AliHFEreducedMCParticle)

//_______________________________________
AliHFEreducedMCParticle::AliHFEreducedMCParticle():
TObject(),
  fSignedPt(0.),
  fP(0.),
  fEta(0.),
  fPhi(0.),
  fPdg(0),
  fMotherPdg(0),
  fSource(5),
  fSignal(kFALSE)
{
  //
  // Default constructor
  //
        memset(fProductionVertex, 0, sizeof(Double_t) *3);      // Initialise production vertex array
}

//_______________________________________
AliHFEreducedMCParticle::AliHFEreducedMCParticle(const AliHFEreducedMCParticle &ref):
  TObject(ref),
  fSignedPt(ref.fSignedPt),
  fP(ref.fP),
  fEta(ref.fEta),
  fPhi(ref.fPhi),
  fPdg(ref.fPdg),
  fMotherPdg(ref.fMotherPdg),
  fSource(ref.fSource),
  fSignal(ref.fSignal)
{
  //
  // Copy constructor
  //
  memcpy(fProductionVertex, ref.fProductionVertex, sizeof(Double_t) *3);      // Port production vertex array
}

//_______________________________________
AliHFEreducedMCParticle &AliHFEreducedMCParticle::operator=(const AliHFEreducedMCParticle &ref){
  // 
  // Assignment operator
  //
  if(&ref != this){
    TObject::operator=(ref);
    fSignedPt = ref.fSignedPt;
    fP= ref.fP;
    fEta = ref.fEta;
    fPhi = ref.fPhi;
    fPdg = ref.fPdg;
    fMotherPdg = ref.fMotherPdg;
    fSource = ref.fSource;
    fSignal = ref.fSignal;
    memcpy(fProductionVertex, ref.fProductionVertex, sizeof(Double_t) *3); 
  }
  return *this;
}

