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
// Class for TOF PID
// Implements the abstract base class AliHFEpidBase
// IsInitialized() does the PID decision
// 
// Authors:
//   Markus Fasel  <M.Fasel@gsi.de>
//   Matus Kalisky <matus.kalisky@cern.ch>  (contact)
//

#include <TMath.h>

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliHFEdetPIDqa.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpidQAmanager.h"


ClassImp(AliHFEpidTOF)

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF():
  AliHFEpidBase()
  , fNsigmaTOF(3)
{
  //
  // Constructor
  //
} 

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const Char_t *name):
  AliHFEpidBase(name)
  , fNsigmaTOF(3)
{
  //
  // Constructor
  //
}

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const AliHFEpidTOF &c):
  AliHFEpidBase("")
  , fNsigmaTOF(3)
{  
  // 
  // Copy operator
  //

  c.Copy(*this);
}
//___________________________________________________________________
AliHFEpidTOF &AliHFEpidTOF::operator=(const AliHFEpidTOF &ref){
  //
  // Assignment operator
  //

  if(this != &ref){
    ref.Copy(*this);
  }

  return *this;
}
//___________________________________________________________________
AliHFEpidTOF::~AliHFEpidTOF(){
  //
  // Destructor
  //
}
//___________________________________________________________________
void AliHFEpidTOF::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidTOF &target = dynamic_cast<AliHFEpidTOF &>(ref);

  target.fNsigmaTOF = fNsigmaTOF;

  AliHFEpidBase::Copy(ref);
}
//___________________________________________________________________
Bool_t AliHFEpidTOF::InitializePID(Int_t /*run*/){
  //
  // InitializePID: TOF experts have to implement code here
  //
  return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidTOF::IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const
{
  //
  // TOF PID based on n-Sigma cut
  // Selects Protons and Kaons via n-sigma cut up to 3 GeV/c
  // In addition histos for n-sigma before (all species) and after (only closest species) are filled
  //
  if(!fkPIDResponse) return 0;
  AliDebug(2, "PID object available");

  const AliVTrack *vtrack = dynamic_cast<const AliVTrack *>(track->GetRecTrack());
  if(!(vtrack && (vtrack->GetStatus() & AliESDtrack::kTOFpid))) return 0;
  AliDebug(2, "Track Has TOF PID");

  if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTOFpid, AliHFEdetPIDqa::kBeforePID);

  // Fill before selection
  Double_t sigEle = fkPIDResponse->NumberOfSigmasTOF(track->GetRecTrack(), AliPID::kElectron);
  AliDebug(2, Form("Number of sigmas in TOF: %f", sigEle));
  Int_t pdg = 0;
  if(TMath::Abs(sigEle) < fNsigmaTOF){
    pdg = 11;
    if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTOFpid, AliHFEdetPIDqa::kAfterPID);
  }
 
  return pdg;
}

