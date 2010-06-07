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
#include "TMath.h"
#include "TText.h"
#include "TLine.h"
#include <TClonesArray.h>

#include <stdio.h>
#include <math.h>

#include "AliPMDUtility.h"
#include "AliAlignObjMatrix.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliLog.h"


ClassImp(AliPMDUtility)

AliPMDUtility::AliPMDUtility():
  fAlObj(GetAlignObj()),
  fPx(0.),
  fPy(0.),
  fPz(0.),
  fTheta(0.),
  fEta(0.),
  fPhi(0.),
  fWriteModule(1)
{
  // Default constructor
  for (Int_t i = 0; i < 4; i++)
    {
      for (Int_t j = 0; j < 3; j++)
	{
	  fSecTr[i][j] = 0.;
	}
    }

}

AliPMDUtility::AliPMDUtility(Float_t px, Float_t py, Float_t pz):
  fAlObj(GetAlignObj()),
  fPx(px),
  fPy(py),
  fPz(pz),
  fTheta(0.),
  fEta(0.),
  fPhi(0.),
  fWriteModule(1)
{
  // Constructor
  for (Int_t i = 0; i < 4; i++)
    {
      for (Int_t j = 0; j < 3; j++)
	{
	  fSecTr[i][j] = 0.;
	}
    }

}
AliPMDUtility::AliPMDUtility(const AliPMDUtility &pmdutil):
  TObject(pmdutil),
  fAlObj(pmdutil.GetAlignObj()),
  fPx(pmdutil.fPx),
  fPy(pmdutil.fPy),
  fPz(pmdutil.fPz),
  fTheta(pmdutil.fTheta),
  fEta(pmdutil.fEta),
  fPhi(pmdutil.fPhi),
  fWriteModule(pmdutil.fWriteModule)
{
  // copy constructor
    for (Int_t i = 0; i < 4; i++)
    {
      for (Int_t j = 0; j < 3; j++)
	{
	  fSecTr[i][j] = pmdutil.fSecTr[i][j];
	}
    }

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
      fWriteModule = pmdutil.fWriteModule;
      for (Int_t i = 0; i < 4; i++)
	{
	  for (Int_t j = 0; j < 3; j++)
	    {
	      fSecTr[i][j] = pmdutil.fSecTr[i][j];
	    }
	}

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
  // Apply the alignment here to the x, y values
  if(ism < 6)
    {
      xpos += fSecTr[0][0];
      ypos += fSecTr[0][1];
    }
  else if(ism >= 6 && ism < 12)
    {
      xpos += fSecTr[1][0];
      ypos += fSecTr[1][1];
    }
  else if(ism >=12 && ism < 18)
    {
      xpos += fSecTr[2][0];
      ypos += fSecTr[2][1];
    }
  else if(ism >= 18 && ism < 24)
    {
      xpos += fSecTr[3][0];
      ypos += fSecTr[3][1];
    }

}
// ---------------------------------------------------------- 
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

  // Apply the alignment here to the x, y values
  if(ism < 6)
    {
      xpos += fSecTr[0][0];
      ypos += fSecTr[0][1];
    }
  else if(ism >= 6 && ism < 12)
    {
      xpos += fSecTr[1][0];
      ypos += fSecTr[1][1];
    }
  else if(ism >=12 && ism < 18)
    {
      xpos += fSecTr[2][0];
      ypos += fSecTr[2][1];
    }
  else if(ism >= 18 && ism < 24)
    {
      xpos += fSecTr[3][0];
      ypos += fSecTr[3][1];
    }

}

// -------------------------------------------------------- //

void AliPMDUtility::RectGeomCellPos(Int_t ism, Float_t xpad,
				    Float_t ypad, Float_t &xpos,
				    Float_t &ypos, Float_t & zpos)
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

  // Apply the alignment here to the x, y, and z values
  if(ism < 6)
    {
      xpos += fSecTr[0][0];
      ypos += fSecTr[0][1];
      zpos += fSecTr[0][2];
    }
  else if(ism >= 6 && ism < 12)
    {
      xpos += fSecTr[1][0];
      ypos += fSecTr[1][1];
      zpos += fSecTr[1][2];
    }
  else if(ism >=12 && ism < 18)
    {
      xpos += fSecTr[2][0];
      ypos += fSecTr[2][1];
      zpos += fSecTr[2][2];
    }
  else if(ism >= 18 && ism < 24)
    {
      xpos += fSecTr[3][0];
      ypos += fSecTr[3][1];
      zpos += fSecTr[3][2];
    }



}
// -------------------------------------------------------- //

void AliPMDUtility::GenerateBoundaryPoints(Int_t ism, Float_t &x1ism, 
					   Float_t &y1ism, Float_t &x2ism,
					   Float_t &y2ism)
{
  // Generate bounding-box.


    Float_t xism = 0, yism = 0;
    Float_t dxism = 0., dyism = 0.;

    const Float_t kRad     = 0.25;
    const Float_t kSqRoot3 = 1.732050808;
    const Float_t kDia     = 0.50;


  const Double_t kXcorner[24] =
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


  const Double_t kYcorner[24] =
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


  if (ism > 23) ism -= 24;


  if (ism < 6)
    {
      xism  = kXcorner[ism] + kRad;
      yism  = kYcorner[ism] + kRad;
      dxism = -kRad*kSqRoot3*48.;
      dyism = -kDia*96. - kRad;
  }
  if (ism >= 6 && ism < 12)
    {
      xism  = kXcorner[ism] - kRad;
      yism  = kYcorner[ism] - kRad;
      dxism = kRad*kSqRoot3*48.;
      dyism = kDia*96. + kRad;
  }
  if (ism >= 12 && ism < 18)
    {
      xism  = kXcorner[ism] + kRad;
      yism  = kYcorner[ism] + kRad;
      dxism = -kRad*kSqRoot3*96.;
      dyism = -kDia*48. - kRad;
  }
  if (ism >= 18 && ism < 24)
    {
      xism  = kXcorner[ism] - kRad;
      yism  = kYcorner[ism] - kRad;
      dxism = kRad*kSqRoot3*96.;
      dyism = kDia*48. + kRad;
  }

  x1ism = xism;
  x2ism = xism + dxism;
  y1ism = yism;
  y2ism = yism + dyism;

}
// ------------------------------------------------------------------- //

void AliPMDUtility::DrawPMDModule(Int_t idet)
{

    Float_t x1ism, x2ism, y1ism, y2ism;
    Float_t deltaX, deltaY;
    
    //TH2F *h2 = new TH2F("h2","Y vs. X",200,-100.,100.,200,-100.,100.);
    //h2->Draw();

    TLine t;
    t.SetLineColor(2);

    TText tt;
    tt.SetTextColor(4);

    Char_t smnumber[10];

    for(Int_t ism=0; ism < 24; ism++)
    {
	GenerateBoundaryPoints(ism, x1ism, y1ism, x2ism, y2ism);
	deltaX = (x2ism - x1ism)/2.;
	deltaY = (y2ism - y1ism)/2.;
	if (fWriteModule == 1)
	{
	  if(idet == 0)
	    {
	      sprintf(smnumber,"%d",ism);
	    }
	  else if (idet == 1)
	    {
	      sprintf(smnumber,"%d",24+ism);
	    }
	    tt.DrawText(x1ism+deltaX,y1ism+deltaY,smnumber);
	}
	t.DrawLine(x1ism, y1ism, x1ism, y2ism);
	t.DrawLine(x1ism, y1ism, x2ism, y1ism);
	t.DrawLine(x2ism, y1ism, x2ism, y2ism);
	t.DrawLine(x1ism, y2ism, x2ism, y2ism);
    }

}

// ------------------------------------------------------------------- //


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
  // Get the alignment stuff here

  AliAlignObjMatrix * aam;
  Double_t tr[3];
  //Double_t secTr[4][3];

  for (Int_t isector=0; isector<4; isector++)
    {
      aam = (AliAlignObjMatrix*)fAlObj->UncheckedAt(isector);
      aam->GetTranslation(tr);
      
      for(Int_t ixyz=0; ixyz < 3; ixyz++)
	{
	  fSecTr[isector][ixyz] = (Float_t) tr[ixyz];
	}
    }
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
void AliPMDUtility::SetWriteModule(Int_t wrmod)
{
    fWriteModule = wrmod;
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
//--------------------------------------------------------------------//
TClonesArray* AliPMDUtility::GetAlignObj() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Align/Data");
  
  if(!entry) AliFatal("Alignment object retrieval failed!");
  
  TClonesArray *alobj = 0;
  if (entry) alobj = (TClonesArray*) entry->GetObject();
  
  if (!alobj)  AliFatal("No alignment data from  database !");
  
  return alobj;
}


