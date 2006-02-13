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

/*
$Log$
Revision 1.1  2006/02/13 16:10:48  arcelli
Add classes for TOF Calibration (C.Zampolli)

author: Chiara Zampolli, zampolli@bo.infn.it
*/  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFcalibESD.h"
#include <Riostream.h>
#include <stdlib.h>


ClassImp(AliTOFcalibESD)

//________________________________________________________________
AliTOFcalibESD::AliTOFcalibESD():
  fToT(0),
  fIntLen(0),
  fTOFtime(0),
  fP(0),
  fSector(-1),
  fPlate(-1),
  fStrip(-1),
  fPadz(-1),
  fPadx(-1),
  fTOFsignalND(0)
{
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i] = 0;
}
//________________________________________________________________

AliTOFcalibESD::AliTOFcalibESD(const AliTOFcalibESD& UnCalib):
  AliESDtrack(UnCalib),
  fToT(UnCalib.fToT),
  fIntLen(UnCalib.fIntLen),
  fTOFtime(UnCalib.fTOFtime),
  fP(UnCalib.fP),
  fCombID(UnCalib.fCombID),
  fSector(UnCalib.fSector),
  fPlate(UnCalib.fPlate),
  fStrip(UnCalib.fStrip),
  fPadz(UnCalib.fPadz),
  fPadx(UnCalib.fPadx),
  fTOFsignalND(UnCalib.fTOFsignalND)
{
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i] = UnCalib.fTrackTime[i];
  for (Int_t i = 0;i<15;i++) fExtCov[i] = UnCalib.fExtCov[i];
}
//________________________________________________________________

AliTOFcalibESD::~AliTOFcalibESD()
{
}
//________________________________________________________________

void AliTOFcalibESD::SetIntegratedTime(const Double_t *tracktime){
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i] = tracktime[i];
}
//________________________________________________________________

void AliTOFcalibESD::CopyFromAliESD(const AliESDtrack* track)
{
  fP = track->GetP();
  fTOFtime = track->GetTOFsignal();
  fIntLen = track->GetIntegratedLength();
  Double_t exptime[10]; 
  track->GetIntegratedTimes(exptime);
  for (Int_t i=0;i<AliPID::kSPECIES;i++) {
    fTrackTime[i] = exptime[i];
  }
  Double_t c[15];
  track->GetExternalCovariance(c);
  for (Int_t i = 0;i<15;i++){
    fExtCov[i] = c[i];
  }
}
//______________________________________________________________________

void AliTOFcalibESD::GetExternalCovariance(Double_t cov[15]) const {
  for (Int_t i=0; i<15; i++) cov[i]=fExtCov[i];
}
//______________________________________________________________________

void AliTOFcalibESD::GetIntegratedTimes(Double_t exp[AliPID::kSPECIES]) const {
  for (Int_t i=0; i<AliPID::kSPECIES; i++) exp[i]=fTrackTime[i];
}
//______________________________________________________________________

Int_t AliTOFcalibESD::Compare(const TObject* uncobj) const{
  Double_t c1[15]; 
  this->GetExternalCovariance(c1);
  Double_t c2[15]; 
  ((AliTOFcalibESD*)uncobj)->GetExternalCovariance(c2);
  if (c1[0]*c1[2] <c2[0]*c2[2]) return -1;
  if (c1[0]*c1[2]>c2[0]*c2[2]) return 1;
  return 0;
}


