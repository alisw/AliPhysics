/***************************************************************************
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
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Utility code for ALICE-PMD                         //
//                                                     //
//-----------------------------------------------------//

#include "Riostream.h"
#include "AliPMDUtility.h"
#include "TMath.h"
#include <stdio.h>
#include <math.h>


ClassImp(AliPMDUtility)

AliPMDUtility::AliPMDUtility():
  fPx(0.),
  fPy(0.),
  fPz(0.),
  fTheta(0.),
  fEta(0.),
  fPhi(0.)
{
  // Default constructor
}

AliPMDUtility::AliPMDUtility(Float_t px, Float_t py, Float_t pz):
  fPx(px),
  fPy(py),
  fPz(pz),
  fTheta(0.),
  fEta(0.),
  fPhi(0.)
{
  // Constructor
}
AliPMDUtility::AliPMDUtility(const AliPMDUtility &pmdutil):
  fPx(pmdutil.fPx),
  fPy(pmdutil.fPy),
  fPz(pmdutil.fPz),
  fTheta(pmdutil.fTheta),
  fEta(pmdutil.fEta),
  fPhi(pmdutil.fPhi)
{
  // copy constructor
}
AliPMDUtility & AliPMDUtility::operator=(const AliPMDUtility &pmdutil)
{
  // assignment operator
  if(this != &pmdutil)
    {
      fPx = pmdutil.fPx;
      fPy = pmdutil.fPy;
      fPz = pmdutil.fPz;
      fTheta = pmdutil.fTheta;
      fEta = pmdutil.fEta;
      fPhi = pmdutil.fPhi;
    }
  return *this;
}
AliPMDUtility::~AliPMDUtility()
{
  // Default destructor
}

void AliPMDUtility::RectGeomCellPos(Int_t ism, Int_t xpad, Int_t ypad, Float_t &xpos, Float_t &ypos)
{
  // This routine finds the cell eta,phi for the new PMD rectangular 
  // geometry in ALICE
  // Authors : Bedanga Mohanty and Dipak Mishra - 29.4.2003
  // modified by B. K. Nandi for change of coordinate sys
  //
  // SMA  ---> Supermodule Type A           ( SM - 0)
  // SMAR ---> Supermodule Type A ROTATED   ( SM - 1)
  // SMB  ---> Supermodule Type B           ( SM - 2)
  // SMBR ---> Supermodule Type B ROTATED   ( SM - 3)
  //
  // ism   : Serial module number from 0 to 23 for each plane

 
  // Corner positions (x,y) of the 24 unit moudles in ALICE PMD

  double xcorner[24] =
    {
      74.8833,  53.0045, 31.1255,    //Type-A
      74.8833,  53.0045, 31.1255,    //Type-A
      -74.8833, -53.0044, -31.1255,  //Type-AR
      -74.8833, -53.0044, -31.1255,  //Type-AR
      8.9165, -33.7471,            //Type-B
      8.9165, -33.7471,            //Type-B
      8.9165, -33.7471,            //Type-B
      -8.9165, 33.7471,            //Type-BR
      -8.9165, 33.7471,            //Type-BR
      -8.9165, 33.7471,            //Type-BR
    };

  
  double ycorner[24] =
    {
      86.225,  86.225,  86.225,      //Type-A
      37.075,  37.075,  37.075,      //Type-A
      -86.225, -86.225, -86.225,     //Type-AR
      -37.075, -37.075, -37.075,     //Type-AR
      86.225,  86.225,               //Type-B
      61.075,  61.075,               //Type-B
      35.925,  35.925,               //Type-B
      -86.225, -86.225,              //Type-BR
      -61.075, -61.075,              //Type-BR
      -35.925, -35.925               //Type-BR
    };

  
  const Float_t kSqroot3      = 1.73205;  // sqrt(3.);
  const Float_t kCellRadius   = 0.25;
  
  //
  //Every even row of cells is shifted and placed
  //in geant so this condition
  //
  Float_t cellRadius = 0.25;
  Float_t shift = 0.0;
  if(xpad%2 == 0)
    {
      shift = -cellRadius/2.0;
    }
  else
    {
      shift = 0.0;
    }


  if(ism < 6)
    {
      ypos = ycorner[ism] - (Float_t) xpad*kCellRadius*2.0 + shift;
      xpos = xcorner[ism] - (Float_t) ypad*kSqroot3*kCellRadius;
    }
  else if(ism >=6 && ism < 12)
    {
      ypos = ycorner[ism] + (Float_t) xpad*kCellRadius*2.0 + shift;
      xpos = xcorner[ism] + (Float_t) ypad*kSqroot3*kCellRadius;
    }
  else if(ism >= 12 && ism < 18)
    {
      ypos = ycorner[ism] - (Float_t) xpad*kCellRadius*2.0 + shift;
      xpos = xcorner[ism] - (Float_t) ypad*kSqroot3*kCellRadius;
    }
  else if(ism >= 18 && ism < 24)
    {
      ypos = ycorner[ism] + (Float_t) xpad*kCellRadius*2.0 + shift;
      xpos = xcorner[ism] + (Float_t) ypad*kSqroot3*kCellRadius;
    }

}

void AliPMDUtility::RectGeomCellPos(Int_t ism, Float_t xpad, Float_t ypad, Float_t &xpos, Float_t &ypos)
{
  // If the xpad and ypad inputs are float, then 0.5 is added to it
  // to find the layer which is shifted.
  // This routine finds the cell eta,phi for the new PMD rectangular 
  // geometry in ALICE
  // Authors : Bedanga Mohanty and Dipak Mishra - 29.4.2003
  // modified by B. K. Nnadi for change of coordinate sys
  //
  // SMA  ---> Supermodule Type A           ( SM - 0)
  // SMAR ---> Supermodule Type A ROTATED   ( SM - 1)
  // SMB  ---> Supermodule Type B           ( SM - 2)
  // SMBR ---> Supermodule Type B ROTATED   ( SM - 3)
  //
  // ism   : Serial Module number from 0 to 23 for each plane

  // Corner positions (x,y) of the 24 unit moudles in ALICE PMD

  double xcorner[24] =
    {
      74.8833,  53.0045, 31.1255,    //Type-A
      74.8833,  53.0045, 31.1255,    //Type-A
      -74.8833, -53.0044, -31.1255,  //Type-AR
      -74.8833, -53.0044, -31.1255,  //Type-AR
      8.9165, -33.7471,            //Type-B
      8.9165, -33.7471,            //Type-B
      8.9165, -33.7471,            //Type-B
      -8.9165, 33.7471,            //Type-BR
      -8.9165, 33.7471,            //Type-BR
      -8.9165, 33.7471,            //Type-BR
    };

  

  double ycorner[24] =
    {
      86.225,  86.225,  86.225,      //Type-A
      37.075,  37.075,  37.075,      //Type-A
      -86.225, -86.225, -86.225,     //Type-AR
      -37.075, -37.075, -37.075,     //Type-AR
      86.225,  86.225,               //Type-B
      61.075,  61.075,               //Type-B
      35.925,  35.925,               //Type-B
      -86.225, -86.225,              //Type-BR
      -61.075, -61.075,              //Type-BR
      -35.925, -35.925               //Type-BR
    };


  const Float_t kSqroot3    = 1.73205;  // sqrt(3.);
  const Float_t kCellRadius = 0.25;
  
  //
  //Every even row of cells is shifted and placed
  //in geant so this condition
  //
  Float_t cellRadius = 0.25;
  Float_t shift = 0.0;
  Int_t iirow = (Int_t) (xpad+0.5);
  if(iirow%2 == 0)
    {
      shift = -cellRadius/2.0;
    }
  else
    {
      shift = 0.0;
    }

  if(ism < 6)
    {
      ypos = ycorner[ism] - xpad*kCellRadius*2.0 + shift;
      xpos = xcorner[ism] - ypad*kSqroot3*kCellRadius;
    }
  else if(ism >=6 && ism < 12)
    {
      ypos = ycorner[ism] + xpad*kCellRadius*2.0 + shift;
      xpos = xcorner[ism] + ypad*kSqroot3*kCellRadius;
    }
  else if(ism >= 12 && ism < 18)
    {
      ypos = ycorner[ism] - xpad*kCellRadius*2.0 + shift;
      xpos = xcorner[ism] - ypad*kSqroot3*kCellRadius;
    }
  else if(ism >= 18 && ism < 24)
    {
      ypos = ycorner[ism] + xpad*kCellRadius*2.0 + shift;
      xpos = xcorner[ism] + ypad*kSqroot3*kCellRadius;
    }

}
void AliPMDUtility::ApplyVertexCorrection(Float_t vertex[], Float_t xpos,
					  Float_t ypos, Float_t zpos)
{
  // Not implemented
  fPx = xpos - vertex[0];
  fPy = ypos - vertex[1];
  fPz = zpos - vertex[2];
}
void AliPMDUtility::ApplyAlignment()
{
  // Not implemented
}

void AliPMDUtility::SetPxPyPz(Float_t px, Float_t py, Float_t pz)
{
  fPx = px;
  fPy = py;
  fPz = pz;
}

void AliPMDUtility::SetXYZ(Float_t xpos, Float_t ypos, Float_t zpos)
{
  fPx = xpos;
  fPy = ypos;
  fPz = zpos;
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

      if(fPx > 0 && fPy > 0) phi = phi1;        // 1st Quadrant
      if(fPx < 0 && fPy > 0) phi = 180 - phi1;  // 2nd Quadrant
      if(fPx < 0 && fPy < 0) phi = 180 + phi1;  // 3rd Quadrant
      if(fPx > 0 && fPy < 0) phi = 360 - phi1;  // 4th Quadrant

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
  
  if(fPx == 0)
    {
      if(fPy>0) phi = 90.;
      if(fPy<0) phi = 270.;
    }
  if(fPx != 0)
    {
      pybypx = fPy/fPx;
      if(pybypx < 0) pybypx = - pybypx;
      phi1 = TMath::ATan(pybypx)*180./3.14159;
      if(fPx > 0 && fPy > 0) phi = phi1;        // 1st Quadrant
      if(fPx < 0 && fPy > 0) phi = 180 - phi1;  // 2nd Quadrant
      if(fPx < 0 && fPy < 0) phi = 180 + phi1;  // 3rd Quadrant
      if(fPx > 0 && fPy < 0) phi = 360 - phi1;  // 4th Quadrant

    }
  phi = phi*3.14159/180.;

  fTheta = theta;
  fEta   = eta;
  fPhi   = phi;
}
void AliPMDUtility::CalculateXY(Float_t eta, Float_t phi, Float_t zpos)
{
  // Not implemented

  //  eta   = -TMath::Log(TMath::Tan(0.5*theta));

  Float_t xpos = 0., ypos = 0.;

  //  Float_t theta = 2.0*TMath::ATan(TMath::Log(-eta));

  fEta = eta;
  fPhi = phi;
  fPx  = xpos;
  fPy  = ypos;
  fPz  = zpos;
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
Float_t AliPMDUtility::GetX() const
{
  return fPx;
}
Float_t AliPMDUtility::GetY() const
{
  return fPy;
}
Float_t AliPMDUtility::GetZ() const
{
  return fPz;
}

