/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD Pid class to store additional pid information
//     Author: Annalisa Mastroserio
//-------------------------------------------------------------------------

#include "AliAODPid.h"
#include "AliESDtrack.h"
#include "AliLog.h"

ClassImp(AliAODPid)


//______________________________________________________________________________
AliAODPid::AliAODPid():
    fITSsignal(0), 
    fTPCsignal(0),
    fTRDnSlices(0),
    fTRDslices(0x0),
    fTOFesdsignal(0),
    fHMPIDsignal(0)
{
  // default constructor
    for(Int_t i=0; i<kSPECIES; i++) fIntTime[i]=0; 
}

//______________________________________________________________________________
AliAODPid::~AliAODPid() 
{
  // destructor
}


//______________________________________________________________________________
AliAODPid::AliAODPid(const AliAODPid& pid) : 
  TObject(pid),
  fITSsignal(pid.fITSsignal), 
  fTPCsignal(pid.fTPCsignal),
  fTRDnSlices(pid.fTRDnSlices),
  fTRDslices(0x0),
  fTOFesdsignal(pid.fTOFesdsignal),
  fHMPIDsignal(pid.fHMPIDsignal)
{
  // Copy constructor
    fTRDslices = new Double32_t[fTRDnSlices];
    for(Int_t i=0; i< fTRDnSlices; i++) fTRDslices[i]=pid.fTRDslices[i];
    for(Int_t i=0; i<kSPECIES; i++) fIntTime[i]=pid.fIntTime[i];
}

//______________________________________________________________________________
AliAODPid& AliAODPid::operator=(const AliAODPid& pid)
{
  // Assignment operator
  if(this!=&pid) {
    // copy stuff
  fITSsignal=pid.fITSsignal; 
  fTPCsignal=pid.fTPCsignal;
  fTRDnSlices=pid.fTRDnSlices;
  for(Int_t i=0; i< fTRDnSlices; i++) fTRDslices[i]=pid.fTRDslices[i];
  fTOFesdsignal=pid.fTOFesdsignal;
  fHMPIDsignal=pid.fHMPIDsignal;
  for(Int_t i=0; i<kSPECIES; i++) fIntTime[i]=pid.fIntTime[i];
  }

  return *this;
}
//_______________________________________________________________________________
void AliAODPid::SetDetectorRawSignals(AliESDtrack *track, Double_t timezero)
{
//
//assignment of the detector signals (AliXXXesdPID inspired)
//
 if(!track){
 AliInfo("no ESD track found. .....exiting");
 return;
 }

 fITSsignal=track->GetITSsignal();
 fTPCsignal=track->GetTPCsignal();
 fTRDnSlices=track->GetNumberOfTRDslices()*kTRDnPlanes;
 track->GetIntegratedTimes(fIntTime);
 fTOFesdsignal=track->GetTOFsignal()-timezero; //TO BE FIXED 
 fHMPIDsignal=track->GetHMPIDsignal();

 fTRDslices=new Double32_t[fTRDnSlices];  
 for(Int_t iSl =0; iSl < track->GetNumberOfTRDslices(); iSl++) {
     for(Int_t iPl =0; iPl<kTRDnPlanes; iPl++) fTRDslices[iPl*track->GetNumberOfTRDslices()+iSl] = track->GetTRDslice(iPl,iSl);
    } 
}
//________________________________________________________________________________
void AliAODPid::GetIntegratedTimes(Double_t timeint[5])
{
 // Returns the array with integrated times for each particle hypothesis
for(Int_t i=0; i<kSPECIES; i++) timeint[i]=fIntTime[i];
}
