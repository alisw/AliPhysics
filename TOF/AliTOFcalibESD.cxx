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
Revision 1.4  2006/04/05 08:35:38  hristov
Coding conventions (S.Arcelli, C.Zampolli)

Revision 1.3  2006/03/28 14:57:13  arcelli
updates to handle new V5 geometry & some re-arrangements

Revision 1.2  2006/02/13 17:22:26  arcelli
just Fixing Log info

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

class AliPID;

ClassImp(AliTOFcalibESD)

//________________________________________________________________
AliTOFcalibESD::AliTOFcalibESD():
  fTOFCalCh(-1),
  fToT(0),
  fIntLen(0),
  fTOFtime(0),
  fMo(0),
  fTOFsignalND(0)
{
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrTime[i] = 0;
}
//________________________________________________________________

AliTOFcalibESD::AliTOFcalibESD(const AliTOFcalibESD& UnCalib):
  AliESDtrack(UnCalib),
  fCombID(UnCalib.fCombID),
  fTOFCalCh(UnCalib.fTOFCalCh),
  fToT(UnCalib.fToT),
  fIntLen(UnCalib.fIntLen),
  fTOFtime(UnCalib.fTOFtime),
  fMo(UnCalib.fMo),
  fTOFsignalND(UnCalib.fTOFsignalND)
{
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrTime[i] = UnCalib.fTrTime[i];
  for (Int_t i = 0;i<15;i++) fExtCov[i] = UnCalib.fExtCov[i];
}
//________________________________________________________________

AliTOFcalibESD::~AliTOFcalibESD()
{
}
//________________________________________________________________

void AliTOFcalibESD::SetIntegratedTime(const Double_t *tracktime){
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrTime[i] = tracktime[i];
}
//________________________________________________________________

void AliTOFcalibESD::CopyFromAliESD(const AliESDtrack* track){
  //copy ESD track info 
  fMo = track->GetP();
  fTOFtime = track->GetTOFsignal();
  fToT = track->GetTOFsignalToT();
  fTOFCalCh = track->GetTOFCalChannel();
  fIntLen = track->GetIntegratedLength();
  Double_t exptime[10]; 
  track->GetIntegratedTimes(exptime);
  for (Int_t i=0;i<AliPID::kSPECIES;i++) {
    fTrTime[i] = exptime[i];
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
  for (Int_t i=0; i<AliPID::kSPECIES; i++) exp[i]=fTrTime[i];
}
//______________________________________________________________________

Int_t AliTOFcalibESD::Compare(const TObject* uncobj) const{
  //To order in momentum
  Double_t c1[15]; 
  this->GetExternalCovariance(c1);
  Double_t c2[15]; 
  ((AliTOFcalibESD*)uncobj)->GetExternalCovariance(c2);
  if (c1[0]*c1[2] <c2[0]*c2[2]) return -1;
  if (c1[0]*c1[2]>c2[0]*c2[2]) return 1;
  return 0;
}


