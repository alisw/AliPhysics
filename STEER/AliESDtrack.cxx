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
#include "../ITS/AliITStrackV2.h"

ClassImp(AliESDtrack)

//_______________________________________________________________________
AliESDtrack::AliESDtrack() : 
fFlags(0), 
fITSncls(0),
fTPCncls(0),
fVertex(kFALSE)
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
Bool_t AliESDtrack::UpdateTrackParams(AliKalmanTrack *t, ULong_t flags) {
  //
  // This function updates track's running parameters 
  //
  switch (flags) {
    
  case kITSin:
  case kITSout: 
  case kITSrefit:
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
  case kTRDin: case kTRDout: case kTRDrefit:
    fTRDncls=t->GetNumberOfClusters();
    fTRDchi2=t->GetChi2();
    fTRDsignal=t->GetPIDsignal();
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
  
  if (flags == kITSin)
   {
     AliITStrackV2* itstrack = dynamic_cast<AliITStrackV2*>(t);
     if (itstrack)
      {
        itstrack->PropagateTo(3.,0.0028,65.19);
        itstrack->PropagateToVertex();
        
        Double_t ralpha=t->GetAlpha();
        Double_t rx;      // X-coordinate of the track reference plane 
        Double_t rp[5];   // external track parameters  
        t->GetExternalParameters(rx,rp);
   
        Double_t phi=TMath::ASin(rp[2]) + ralpha;
        Double_t pt=1./TMath::Abs(rp[4]);
        Double_t r=TMath::Sqrt(rx*rx + rp[0]*rp[0]);
        
        fVertexX=r*TMath::Cos(phi); 
        fVertexY=r*TMath::Sin(phi); 
        fVertexZ=rp[1]; 
        
        fVertexPx = pt*TMath::Cos(phi); 
        fVertexPy = pt*TMath::Sin(phi); 
        fVertexPz = pt*rp[3]; 
        fVertex = kTRUE;
      }
   }
  
  return kTRUE;
}

//_______________________________________________________________________
void AliESDtrack::GetExternalParametersAt(Double_t x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters
  // at the plane x
  //---------------------------------------------------------------------
  Double_t dx=x-fRx;
  Double_t c=fRp[4]/AliKalmanTrack::GetConvConst();
  Double_t f1=fRp[2], f2=f1 + c*dx;
  Double_t r1=sqrt(1.- f1*f1), r2=sqrt(1.- f2*f2);    

  p[0]=fRp[0]+dx*(f1+f2)/(r1+r2);
  p[1]=fRp[1]+dx*(f1+f2)/(f1*r2 + f2*r1)*fRp[3];
  p[2]=fRp[2]+dx*c;
  p[3]=fRp[3];
  p[4]=fRp[4];
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
void AliESDtrack::SetITSpid(const Double_t *p) {  
  for (Int_t i=0; i<kSPECIES; i++) fITSr[i]=p[i];
  SetStatus(AliESDtrack::kITSpid);
}

//_______________________________________________________________________
void AliESDtrack::GetITSpid(Double_t *p) const {
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fITSr[i];
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
Int_t AliESDtrack::GetTPCclusters(Int_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  for (Int_t i=0; i<180; i++) idx[i]=fTPCindex[i];  // MI I prefer some constant
  return fTPCncls;
}

//_______________________________________________________________________
void AliESDtrack::SetTPCpid(const Double_t *p) {  
  for (Int_t i=0; i<kSPECIES; i++) fTPCr[i]=p[i];
  SetStatus(AliESDtrack::kTPCpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTPCpid(Double_t *p) const {
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fTPCr[i];
}

//_______________________________________________________________________
void AliESDtrack::SetTRDpid(const Double_t *p) {  
  for (Int_t i=0; i<kSPECIES; i++) fTRDr[i]=p[i];
  SetStatus(AliESDtrack::kTRDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTRDpid(Double_t *p) const {
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fTRDr[i];
}

//_______________________________________________________________________
void    AliESDtrack::SetTRDpid(Int_t iSpecies, Float_t p)
{
  fTRDr[iSpecies] = p;
}

Float_t AliESDtrack::GetTRDpid(Int_t iSpecies) const
{
  return fTRDr[iSpecies];
}

//_______________________________________________________________________
void AliESDtrack::SetTOFpid(const Double_t *p) {  
  for (Int_t i=0; i<kSPECIES; i++) fTOFr[i]=p[i];
  SetStatus(AliESDtrack::kTOFpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTOFpid(Double_t *p) const {
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fTOFr[i];
}

//_______________________________________________________________________
void AliESDtrack::SetESDpid(const Double_t *p) {  
  for (Int_t i=0; i<kSPECIES; i++) fR[i]=p[i];
  SetStatus(AliESDtrack::kESDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetESDpid(Double_t *p) const {
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fR[i];
}

void AliESDtrack::GetVertexXYZ(Double_t& x,Double_t& y, Double_t&z) const
{
//returns track position in DCA to vertex  
  x = fVertexX;
  y = fVertexY;
  z = fVertexZ;
}
void AliESDtrack::GetVertexPxPyPz(Double_t& px,Double_t& py, Double_t& pz) const
{
//returns track momentum in DCA to vertex  
  px = fVertexPx;
  py = fVertexPy;
  pz = fVertexPz;
}
