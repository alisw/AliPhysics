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
// ITS PID class
// checks ITS PID based on ITS dE/dx truncated mean
//
// Authors: Matus Kalisky <matus.kalisky@cern.ch>
//          Markus Fasel <M.Fasel@gsi.de>
//
#include <TClass.h>
#include <TH2F.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>

//#include "AliAODTrack.h"
//#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliVParticle.h"
#include "AliPIDResponse.h"

#include "AliHFEdetPIDqa.h"
#include "AliHFEpidITS.h"
#include "AliHFEpidQAmanager.h"

ClassImp(AliHFEpidITS)

//___________________________________________________________________
AliHFEpidITS::AliHFEpidITS():
AliHFEpidBase()
, fNsigmaITSlow(-3)
, fNsigmaITShigh(3)
, fMeanShift(0)
{
    //
    // Constructor
    //
    
} 

//___________________________________________________________________
AliHFEpidITS::AliHFEpidITS(const Char_t *name):
AliHFEpidBase(name)
, fNsigmaITSlow(-3)
, fNsigmaITShigh(3)
, fMeanShift(0)
{
    //
    // Default constructor
    //
}

//___________________________________________________________________
AliHFEpidITS::AliHFEpidITS(const AliHFEpidITS &ref):
AliHFEpidBase("")
, fNsigmaITSlow(ref.fNsigmaITSlow)
, fNsigmaITShigh(ref.fNsigmaITShigh)
, fMeanShift(ref.fMeanShift)
{
    //
    // Copy constructor
    //
    ref.Copy(*this);
}

//___________________________________________________________________
AliHFEpidITS &AliHFEpidITS::operator=(const AliHFEpidITS &ref){
    //
    // Assignment operator
    //
    if(this != &ref) ref.Copy(*this);
    return *this;
}

//___________________________________________________________________
AliHFEpidITS::~AliHFEpidITS(){
    //
    // Destructor
    //
}

//___________________________________________________________________
void AliHFEpidITS::Copy(TObject &ref) const {
    //
    // Copy function
    // Provides a deep copy
    //
    AliHFEpidITS &target = dynamic_cast<AliHFEpidITS &>(ref);
    target.fNsigmaITSlow = fNsigmaITSlow;
    target.fNsigmaITShigh = fNsigmaITShigh;
    target.fMeanShift = fMeanShift;
    AliHFEpidBase::Copy(ref);
}

//___________________________________________________________________
Bool_t AliHFEpidITS::InitializePID(Int_t /*run*/){
    //
    // ITS PID initialization
    //
    return kTRUE;
}


//___________________________________________________________________
Int_t AliHFEpidITS::IsSelected(const AliHFEpidObject* track, AliHFEpidQAmanager* pidqa) const {
    //
    // Does PID decision for ITS
    //
    if(!fkPIDResponse) return 0;
    AliDebug(2, "PID object available");
    
    const AliVTrack *vtrack = dynamic_cast<const AliVTrack *>(track->GetRecTrack());
    if(!vtrack) return 0;
    
    if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kITSpid, AliHFEdetPIDqa::kBeforePID);
    
    // Fill before selection
    Int_t pdg = 0;
    Double_t sigEle = GetITSNsigmaCorrected(vtrack);
    AliDebug(2, Form("Number of sigmas in ITS: %f", sigEle));
    if(sigEle > fNsigmaITSlow && sigEle < fNsigmaITShigh) pdg = 11;
  //  if(TMath::Abs(sigEle) < fNsigmaITS) pdg = 11;
    if(pdg == 11 && pidqa) pidqa->ProcessTrack(track, AliHFEpid::kITSpid, AliHFEdetPIDqa::kAfterPID);
    return pdg;
    //  return 11;  // @TODO: Implement ITS PID decision
}

//___________________________________________________________________
Double_t AliHFEpidITS::GetITSNsigmaCorrected(const AliVTrack *track) const {
    //
    // Get the ITS number of sigmas corrected for a possible shift of the mean dE/dx
    //
    return fkPIDResponse->NumberOfSigmasITS(track, AliPID::kElectron) - fMeanShift;
}
//___________________________________________________________________
void AliHFEpidITS::SetITSnSigma(Float_t nSigmalow, Float_t nSigmahigh) {
    //
    //Set nSigma cut
    //
    fNsigmaITSlow = nSigmalow;
    fNsigmaITShigh = nSigmahigh;
}
//___________________________________________________________________





