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

//-----------------------------------------------------------------
//           Implementation of the ESD track class
//   This is the class to deal with during the phisical analysis of data
//
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "TMath.h"

#include "AliESDtrack.h"
#include "AliKalmanTrack.h"

ClassImp(AliESDtrack)

//_______________________________________________________________________
AliESDtrack::AliESDtrack() : 
fFlags(0), 
fITSncls(0),
fTPCncls(0)
{
  //
  // The default ESD constructor 
  //
  for (Int_t i=0; i<kSPECIES; i++) fR[i]=0.;
}

//_______________________________________________________________________
Float_t AliESDtrack::GetMass() const {
  Float_t max=0.;
  Int_t k=-1;
  for (Int_t i=0; i<kSPECIES; i++) {
    if (fR[i]>max) {k=i; max=fR[i];}
  }
  if (k==0) return 0.00051;
  if (k==1) return 0.10566;
  if (k==2||k==-1) return 0.13957;
  if (k==3) return 0.49368;
  if (k==4) return 0.93827;
  Warning("GetMass()","Undefined mass !");
  return 0.13957;
}

//_______________________________________________________________________
Bool_t AliESDtrack::UpdateTrackParams(const AliKalmanTrack *t, ULong_t flags) {
  //
  // This function updates track's running parameters 
  //
  switch (flags) {
  case kITSin: case kITSout: case kITSrefit:
    fITSncls=t->GetNumberOfClusters();
    fITSchi2=t->GetChi2();
    for (Int_t i=0;i<fITSncls;i++) fITSindex[i]=t->GetClusterIndex(i);
    fITSsignal=t->GetPIDsignal();
    break;
  case kTPCin: case kTPCout: case kTPCrefit:
    fTPCncls=t->GetNumberOfClusters();
    fTPCchi2=t->GetChi2();
    for (Int_t i=0;i<fTPCncls;i++) fTPCindex[i]=t->GetClusterIndex(i);
    fTPCsignal=t->GetPIDsignal();
    {Double_t mass=t->GetMass();    // preliminary mass setting 
    if (mass>0.5) fR[4]=1.;         //        used by
    else if (mass<0.4) fR[2]=1.;    // the ITS reconstruction
    else fR[3]=1.;}                 //
    break;
  default: 
    Error("UpdateTrackParams()","Wrong flag !\n");
    return kFALSE;
  }

  SetStatus(flags);
  fLabel=t->GetLabel();

  if (t->IsStartedTimeIntegral()) {
    SetStatus(kTIME);
    Double_t times[10];t->GetIntegratedTimes(times); SetIntegratedTimes(times);
    SetIntegratedLength(t->GetIntegratedLength());
  }

  fRalpha=t->GetAlpha();
  t->GetExternalParameters(fRx,fRp);
  t->GetExternalCovariance(fRc);
  return kTRUE;
}

//_______________________________________________________________________
void AliESDtrack::GetExternalParameters(Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters
  //---------------------------------------------------------------------
  x=fRx;
  for (Int_t i=0; i<5; i++) p[i]=fRp[i];
}

Double_t AliESDtrack::GetP() const {
  //---------------------------------------------------------------------
  // This function returns the track momentum
  //---------------------------------------------------------------------
  Double_t lam=TMath::ATan(fRp[3]);
  Double_t pt=1./TMath::Abs(fRp[4]);
  return pt/TMath::Cos(lam);
}

void AliESDtrack::GetPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  //---------------------------------------------------------------------
  Double_t phi=TMath::ASin(fRp[2]) + fRalpha;
  Double_t pt=1./TMath::Abs(fRp[4]);
  p[0]=pt*TMath::Cos(phi); p[1]=pt*TMath::Sin(phi); p[2]=pt*fRp[3]; 
}

void AliESDtrack::GetXYZ(Double_t *xyz) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  //---------------------------------------------------------------------
  Double_t phi=TMath::ASin(fRp[2]) + fRalpha;
  Double_t r=TMath::Sqrt(fRx*fRx + fRp[0]*fRp[0]);
  xyz[0]=r*TMath::Cos(phi); xyz[1]=r*TMath::Sin(phi); xyz[2]=fRp[1]; 
}

//_______________________________________________________________________
void AliESDtrack::GetExternalCovariance(Double_t c[15]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the cov. matrix
  //---------------------------------------------------------------------
  for (Int_t i=0; i<15; i++) c[i]=fRc[i];
}

//_______________________________________________________________________
void AliESDtrack::GetIntegratedTimes(Double_t *times) const {
  for (Int_t i=0; i<kSPECIES; i++) times[i]=fTrackTime[i];
}

//_______________________________________________________________________
void AliESDtrack::SetIntegratedTimes(const Double_t *times) {
  for (Int_t i=0; i<kSPECIES; i++) fTrackTime[i]=times[i];
}

//_______________________________________________________________________
Int_t AliESDtrack::GetITSclusters(UInt_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  for (Int_t i=0; i<fITSncls; i++) idx[i]=fITSindex[i];
  return fITSncls;
}

//_______________________________________________________________________
Int_t AliESDtrack::GetTPCclusters(UInt_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  for (Int_t i=0; i<fTPCncls; i++) idx[i]=fTPCindex[i];
  return fTPCncls;
}
