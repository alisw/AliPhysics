//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Utility code for ALICE-PMD                         //
//                                                     //
//-----------------------------------------------------//

#include "AliPMDUtility.h"
#include "TMath.h"
#include <stdio.h>

ClassImp(AliPMDUtility)

AliPMDUtility::AliPMDUtility()
{
  fPx    = 0.;
  fPy    = 0.;
  fPz    = 0.;
  fTheta = 0.;
  fEta   = 0.;
  fPhi   = 0.;
}

AliPMDUtility::AliPMDUtility(Float_t Px, Float_t Py, Float_t Pz)
{
  fPx = Px;
  fPy = Py;
  fPz = Pz;
  fTheta = 0.;
  fEta   = 0.;
  fPhi   = 0.;
}

AliPMDUtility::~AliPMDUtility()
{

}

void AliPMDUtility::SetPxPyPz(Float_t Px, Float_t Py, Float_t Pz)
{
  fPx = Px;
  fPy = Py;
  fPz = Pz;
}

void AliPMDUtility::SetXYZ(Float_t xPos, Float_t yPos, Float_t zPos)
{
  fPx = xPos;
  fPy = yPos;
  fPz = zPos;
}
void AliPMDUtility::CalculateEta()
{
  Float_t rpxpy, theta, eta;

  rpxpy  = TMath::Sqrt(fPx*fPx + fPy*fPy);
  theta  = TMath::ATan2(rpxpy,fPz);
  eta    = -TMath::Log(TMath::Tan(0.5*theta));
  fTheta = theta;
  fEta   = eta;
}
void AliPMDUtility::CalculatePhi()
{
  Float_t pybypx, phi = 0., phi1;

  if(fPx==0)
    {
      if(fPy>0) phi = 90.;
      if(fPy<0) phi = 270.;
    }
  if(fPx != 0)
    {
      pybypx = fPy/fPx;
      if(pybypx < 0) pybypx = - pybypx;
      phi1 = TMath::ATan(pybypx)*180./3.14159;
      if(fPx < 0 && fPy > 0) phi = 180 - phi1;
      if(fPx < 0 && fPy < 0) phi = 180 + phi1;
      if(fPx > 0 && fPy < 0) phi = 360 - phi1;
      if(fPx > 0 && fPy > 0) phi = phi1;
    }
  phi = phi*3.14159/180.;

  fPhi = phi;

}
void AliPMDUtility::CalculateEtaPhi()
{
  Float_t rpxpy, theta, eta;
  Float_t pybypx, phi = 0., phi1;

  rpxpy = TMath::Sqrt(fPx*fPx + fPy*fPy);
  theta = TMath::ATan2(rpxpy,fPz);
  eta   = -TMath::Log(TMath::Tan(0.5*theta));
  
  if(fPx==0)
    {
      if(fPy>0) phi = 90.;
      if(fPy<0) phi = 270.;
    }
  if(fPx != 0)
    {
      pybypx = fPy/fPx;
      if(pybypx < 0) pybypx = - pybypx;
      phi1 = TMath::ATan(pybypx)*180./3.14159;
      if(fPx < 0 && fPy > 0) phi = 180 - phi1;
      if(fPx < 0 && fPy < 0) phi = 180 + phi1;
      if(fPx > 0 && fPy < 0) phi = 360 - phi1;
      if(fPx > 0 && fPy > 0) phi = phi1;
    }
  phi = phi*3.14159/180.;

  fTheta = theta;
  fEta   = eta;
  fPhi   = phi;
}
Float_t AliPMDUtility::GetTheta() const
{
  return fTheta;
}
Float_t AliPMDUtility::GetEta() const
{
  return fEta;
}
Float_t AliPMDUtility::GetPhi() const
{
  return fPhi;
}

