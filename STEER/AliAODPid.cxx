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
  delete [] fTRDslices;
  fTRDslices = 0;
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
    
    if(pid.fTRDnSlices<=0||(fTRDnSlices!=pid.fTRDnSlices)){
      // only delete if number changed or is 0
      delete [] fTRDslices;
      fTRDslices = 0;
      if(pid.fTRDnSlices>0) fTRDslices = new Double32_t[fTRDnSlices];
    }
    fTRDnSlices=pid.fTRDnSlices;
    
    for(Int_t i=0; i< fTRDnSlices; i++) fTRDslices[i]=pid.fTRDslices[i];
    fTOFesdsignal=pid.fTOFesdsignal;
    fHMPIDsignal=pid.fHMPIDsignal;
    for(Int_t i=0; i<kSPECIES; i++) fIntTime[i]=pid.fIntTime[i];
  }

  return *this;
}
//_______________________________________________________________________________
void AliAODPid::GetIntegratedTimes(Double_t timeint[kSPECIES])
{
 // Returns the array with integrated times for each particle hypothesis
for(Int_t i=0; i<kSPECIES; i++) timeint[i]=fIntTime[i];
}
//_______________________________________________________________________________
void AliAODPid::SetIntegratedTimes(Double_t timeint[kSPECIES])
{
 // Returns the array with integrated times for each particle hypothesis
for(Int_t i=0; i<kSPECIES; i++) fIntTime[i]=timeint[i];
}
