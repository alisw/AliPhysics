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

#include "AliAODPid.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliTPCdEdxInfo.h"

ClassImp(AliAODPid)


//______________________________________________________________________________
AliAODPid::AliAODPid():
    fITSsignal(0), 
    fTPCsignal(0),
    fTPCsignalN(0),
    fTPCmomentum(0),
    fTPCTgl(0),
    fTRDnSlices(0),
    fTRDntls(0),
    fTRDslices(0x0),
    fTRDsignal(0),
    fTRDChi2(0x0),
    fTOFesdsignal(0),
    fTPCdEdxInfo(0)
{
  // default constructor
    for(Int_t i=0; i<AliPID::kSPECIES; i++) fIntTime[i]   = 0; 
    for(Int_t i=0; i<5; i++) fTOFpidResolution[i] = 0.;
    for(Int_t i=0; i<6; i++) {
      fTRDmomentum[i]      = 0.;
      fTRDncls[i]          = 0;
    }
    for(Int_t i=0; i<4; i++) fITSdEdxSamples[i]   = 0.;
}

//______________________________________________________________________________
AliAODPid::~AliAODPid() 
{
  delete [] fTRDslices;
  fTRDslices = 0;
  delete fTPCdEdxInfo;
  // destructor
}


//______________________________________________________________________________
AliAODPid::AliAODPid(const AliAODPid& pid) : 
  TObject(pid),
  fITSsignal(pid.fITSsignal), 
  fTPCsignal(pid.fTPCsignal),
  fTPCsignalN(pid.fTPCsignalN),
  fTPCmomentum(pid.fTPCmomentum),
  fTPCTgl(pid.fTPCTgl),
  fTRDnSlices(pid.fTRDnSlices),
  fTRDntls(pid.fTRDntls),
  fTRDslices(0x0),
  fTRDsignal(pid.fTRDsignal),
  fTRDChi2(pid.fTRDChi2),
  fTOFesdsignal(pid.fTOFesdsignal),
  fTPCdEdxInfo(0x0)
{
  /// Copy constructor

  SetTRDslices(fTRDnSlices, pid.fTRDslices);
    for(Int_t i=0; i<AliPID::kSPECIES; i++) fIntTime[i]=pid.fIntTime[i];

    for(Int_t i=0; i<6; i++){ 
      fTRDmomentum[i]=pid.fTRDmomentum[i];
      fTRDncls[i] = 0;
    }

    for(Int_t i=0; i<5; i++) fTOFpidResolution[i]=pid.fTOFpidResolution[i];

    for(Int_t i=0; i<4; i++) fITSdEdxSamples[i]=pid.fITSdEdxSamples[i];

    if (pid.fTPCdEdxInfo) fTPCdEdxInfo=new AliTPCdEdxInfo(*pid.fTPCdEdxInfo);
}

//______________________________________________________________________________
AliAODPid& AliAODPid::operator=(const AliAODPid& pid)
{
  /// Assignment operator

  if(this!=&pid) {
    // copy stuff
    TObject::operator=(pid);

    fITSsignal   = pid.fITSsignal; 
    for (Int_t i = 0; i < 4; i++) fITSdEdxSamples[i]=pid.fITSdEdxSamples[i];
    fTPCsignal   = pid.fTPCsignal;
    fTPCsignalN  = pid.fTPCsignalN;
    fTPCmomentum = pid.fTPCmomentum;
    fTPCTgl      = pid.fTPCTgl;

    fTRDsignal = pid.fTRDsignal;
    if(fTRDnSlices != pid.fTRDnSlices) {
      // only delete if number changed or is 0
      delete [] fTRDslices;
      fTRDslices = 0;
      fTRDnSlices = pid.fTRDnSlices;
      if(pid.fTRDnSlices > 0) fTRDslices = new Double32_t[fTRDnSlices];
    }

    if (fTRDslices && pid.fTRDslices)
      memcpy(fTRDslices, pid.fTRDslices, fTRDnSlices*sizeof(Double32_t));

    fTRDntls = pid.fTRDntls;
    for(Int_t i = 0; i < 6; i++){ 
	fTRDmomentum[i] = pid.fTRDmomentum[i];
	fTRDncls[i]     = pid.fTRDncls[i];
    }

    fTRDChi2 = pid.fTRDChi2;

    fTOFesdsignal=pid.fTOFesdsignal;
    for (Int_t i = 0; i < 5; i++) fTOFpidResolution[i]=pid.fTOFpidResolution[i];
    for (Int_t i = 0; i < 5; i++) fIntTime[i]=pid.fIntTime[i];
    
     SetTPCdEdxInfo(pid.fTPCdEdxInfo);
  }

  return *this;
}
//_______________________________________________________________________________
void AliAODPid::GetIntegratedTimes(Double_t *timeint, Int_t nspec) const
{
  /// Returns the array with integrated times for each particle hypothesis

  for(Int_t i=0; i<AliPID::kSPECIES; i++) timeint[i]=fIntTime[i];
  //Note: at the moment only kSPECIES entries are available
  if (nspec>AliPID::kSPECIES) for (int i=AliPID::kSPECIES;i<AliPID::kSPECIESC;i++) timeint[i]=0;
  //
}
//_______________________________________________________________________________
void AliAODPid::SetIntegratedTimes(Double_t timeint[AliPID::kSPECIES])
{
 /// Returns the array with integrated times for each particle hypothesis

 for(Int_t i=0; i<AliPID::kSPECIES; i++) fIntTime[i]=timeint[i];
}
//______________________________________________________________________________
void AliAODPid::SetTOFpidResolution(Double_t tofPIDres[5])
{
  for (Int_t i=0; i<5; i++) fTOFpidResolution[i]=tofPIDres[i];

}
//______________________________________________________________________________
void AliAODPid::GetTOFpidResolution(Double_t tofRes[5]) const
{
  for (Int_t i=0; i<5; i++) tofRes[i]=fTOFpidResolution[i];
}

//______________________________________________________________________________
void AliAODPid::SetITSdEdxSamples(const Double_t s[4])
{
  /// Set the 4 values of dE/dx from individual ITS layers that are read from ESD

  for (Int_t i=0; i<4; i++) fITSdEdxSamples[i]=s[i];
}
//______________________________________________________________________________
void AliAODPid::GetITSdEdxSamples(Double_t s[4]) const
{
  /// Get the 4 values of dE/dx from individual ITS layers that are read from ESD

  for (Int_t i=0; i<4; i++) s[i]=fITSdEdxSamples[i];
}
//______________________________________________________________________________
void AliAODPid::SetTPCdEdxInfo(AliTPCdEdxInfo * dEdxInfo)
{
  /// Set TPC dEdx info

  if (dEdxInfo==0x0){
    delete fTPCdEdxInfo;
    fTPCdEdxInfo=0x0;
    return;
  }
  if (!fTPCdEdxInfo) fTPCdEdxInfo=new AliTPCdEdxInfo;
  (*fTPCdEdxInfo)=(*dEdxInfo);
}

