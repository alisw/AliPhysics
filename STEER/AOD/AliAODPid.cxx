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
#include "AliPID.h"

ClassImp(AliAODPid)


//______________________________________________________________________________
AliAODPid::AliAODPid():
    fITSsignal(0), 
    fTPCsignal(0),
    fTPCsignalN(0),
    fTPCmomentum(0),
    fTRDnSlices(0),
    fTRDntls(0),
    fTRDslices(0x0),
    fTOFesdsignal(0),
    fHMPIDsignal(0)
{
  // default constructor
    for(Int_t i=0; i<kSPECIES; i++) fIntTime[i]   = 0; 
    for(Int_t i=0; i<5; i++) fHMPIDprobs[i] = 0.;
    for(Int_t i=0; i<3; i++) fEMCALPosition[i]    = 0.;
    for(Int_t i=0; i<5; i++) fTOFpidResolution[i] = 0.;
    for(Int_t i=0; i<6; i++) {
      fTRDmomentum[i]      = 0.;
      fTRDncls[i]          = 0;
    }
    for(Int_t i=0; i<3; i++) fEMCALMomentum[i]    = 0.;
    for(Int_t i=0; i<4; i++) fITSdEdxSamples[i]   = 0.;
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
  fTPCsignalN(pid.fTPCsignalN),
  fTPCmomentum(pid.fTPCmomentum),
  fTRDnSlices(pid.fTRDnSlices),
  fTRDntls(pid.fTRDntls),
  fTRDslices(0x0),
  fTOFesdsignal(pid.fTOFesdsignal),
  fHMPIDsignal(pid.fHMPIDsignal)
{
  // Copy constructor
  SetTRDsignal(fTRDnSlices, pid.fTRDslices);
    for(Int_t i=0; i<kSPECIES; i++) fIntTime[i]=pid.fIntTime[i];
    for(Int_t i=0; i<5; i++) fHMPIDprobs[i] = pid.fHMPIDprobs[i];
    for(Int_t i=0; i<3; i++) {
      fEMCALPosition[i]=pid.fEMCALPosition[i];
      fEMCALMomentum[i]=pid.fEMCALMomentum[i];
    }
    for(Int_t i=0; i<6; i++){ 
      fTRDmomentum[i]=pid.fTRDmomentum[i];
      fTRDncls[i] = 0;
    }

    for(Int_t i=0; i<5; i++) fTOFpidResolution[i]=pid.fTOFpidResolution[i];

    for(Int_t i=0; i<4; i++) fITSdEdxSamples[i]=pid.fITSdEdxSamples[i];
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
    for(Int_t i=0; i<5; i++) fHMPIDprobs[i] = pid.fHMPIDprobs[i];
    for(Int_t i=0; i<6; i++){ 
      fTRDmomentum[i]=pid.fTRDmomentum[i];
      fTRDncls[i] = pid.fTRDncls[i];
    }
    for(Int_t i=0; i<3; i++) {
      fEMCALPosition[i]=pid.fEMCALPosition[i];
      fEMCALMomentum[i]=pid.fEMCALMomentum[i];
    }
    for (Int_t i=0; i<5; i++) fTOFpidResolution[i]=pid.fTOFpidResolution[i];
    for (Int_t i=0; i<4; i++) fITSdEdxSamples[i]=pid.fITSdEdxSamples[i];
  }

  return *this;
}
//_______________________________________________________________________________
void AliAODPid::GetIntegratedTimes(Double_t timeint[kSPECIES]) const
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
//_______________________________________________________________________________
void AliAODPid::GetEMCALPosition(Double_t emcalpos[3]) const
{
 // Returns the array with extrapolated track position at the EMCAL surface
  for(Int_t i=0; i<3; i++) emcalpos[i]=fEMCALPosition[i];
}
//_______________________________________________________________________________
void AliAODPid::SetEMCALPosition(Double_t emcpos[3])
{
 // Sets the array with extrapolated track position at the EMCAL surface
  for(Int_t i=0; i<3; i++) fEMCALPosition[i]=emcpos[i];
}
//_______________________________________________________________________________
void AliAODPid::GetEMCALMomentum(Double_t emcalmom[3]) const
{
 // Returns the array with extrapolated track momentum at the EMCAL surface
  for(Int_t i=0; i<3; i++) emcalmom[i]=fEMCALMomentum[i];
}
//_______________________________________________________________________________
void AliAODPid::SetEMCALMomentum(Double_t emcmom[3])
{
 // Sets the array with extrapolated track momentum at the EMCAL surface
  for(Int_t i=0; i<3; i++) fEMCALMomentum[i]=emcmom[i];
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
void AliAODPid::SetHMPIDprobs(Double_t hmpPid[5]) 
{
  //
  // Set the HMPID PID probablities that are read from ESD
  //  
  for(Int_t i = 0; i < 5; i++ ) fHMPIDprobs[i] =  hmpPid[i];
}
//______________________________________________________________________________
void AliAODPid::GetHMPIDprobs(Double_t *p) const
{
  //
  // Set the HMPID PID probablities that are read from ESD
  //  
  for(Int_t i = 0; i < AliPID::kSPECIES; i++ ) p[i] =  fHMPIDprobs[i];
}
//______________________________________________________________________________
void AliAODPid::SetITSdEdxSamples(const Double_t s[4])
{
  //
  // Set the 4 values of dE/dx from individual ITS layers that are read from ESD
  //  
  for (Int_t i=0; i<4; i++) fITSdEdxSamples[i]=s[i];
}
//______________________________________________________________________________
void AliAODPid::GetITSdEdxSamples(Double_t s[4]) const
{
  //
  // Get the 4 values of dE/dx from individual ITS layers that are read from ESD
  //  
  for (Int_t i=0; i<4; i++) s[i]=fITSdEdxSamples[i];
}
