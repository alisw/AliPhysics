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
Revision 1.6  2000/10/02 21:28:06  fca
Removal of useless dependecies via forward declarations

Revision 1.5  2000/06/09 20:37:20  morsch
All coding rule violations except RS3 corrected

Revision 1.4  1999/11/03 17:43:20  fca
New version from G.Martinez & A.Morsch

Revision 1.3  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log

*/

#include "AliGenScan.h"
#include "AliRun.h"

 ClassImp(AliGenScan)
    
 AliGenScan::AliGenScan()
	 :AliGenerator(-1)
{
// Constructor
    fXmin=0;
    fXmax=0;
    fNx=1;
    fYmin=0;
    fYmax=0;
    fNy=1;
    fZmin=0;
    fZmax=0;
    fNz=1;
//
//  Read all particles
    fNpart=-1;
}

AliGenScan::AliGenScan(Int_t npart)
    :AliGenerator(npart)
{
// Constructor
    fXmin=0;
    fXmax=0;
    fNx=1;
    fYmin=0;
    fYmax=0;
    fNy=1;
    fZmin=0;
    fZmax=0;
    fNz=1;
}

//____________________________________________________________
AliGenScan::~AliGenScan()
{
// Destructor
}

void AliGenScan::SetRange(Int_t nx, Float_t xmin, Float_t xmax,
		     Int_t ny, Float_t ymin, Float_t ymax,
		     Int_t nz, Float_t zmin, Float_t zmax)
{
// Define the grid
    fXmin=xmin;
    fXmax=xmax;
    fNx=nx;
    fYmin=ymin;
    fYmax=ymax;
    fNy=ny;
    fZmin=zmin;
    fZmax=zmax;
    fNz=nz;
}

//____________________________________________________________
void AliGenScan::Generate()
{
  //
  // Generate one trigger
  //
  
  Float_t polar[3]= {0,0,0};
  //
  Float_t origin[3];
  Float_t p[3];
  Int_t nt;
  Float_t pmom, theta, phi;
  //
  Float_t random[6];
  Float_t dx,dy,dz;
  
  //
  if (fNy > 0) {
      dx=(fXmax-fXmin)/fNx;
  } else {
      dx=1e10;
  }

  if (fNy > 0) {
      dy=(fYmax-fYmin)/fNy;
  } else {
      dy=1e10;
  }

  if (fNz > 0) {
      dz=(fZmax-fZmin)/fNz;
  } else {
      dz=1e10;
  }
  for (Int_t ix=0; ix<fNx; ix++) {
      for (Int_t iy=0; iy<fNy; iy++) {
	  for (Int_t iz=0; iz<fNz; iz++){
	      Rndm(random,6);
	      origin[0]=fXmin+ix*dx+2*(random[0]-0.5)*fOsigma[0];
	      origin[1]=fYmin+iy*dy+2*(random[1]-0.5)*fOsigma[1];
	      origin[2]=fZmin+iz*dz+2*(random[2]-0.5)*fOsigma[2];	     
	      pmom=fPMin+random[3]*(fPMax-fPMin);
	      theta=fThetaMin+random[4]*(fThetaMax-fThetaMin);
	      phi=fPhiMin+random[5]*(fPhiMax-fPhiMin);
	      p[0] = pmom*TMath::Cos(phi)*TMath::Sin(theta);
	      p[1] = pmom*TMath::Sin(phi)*TMath::Sin(theta);
	      p[2] = pmom*TMath::Cos(theta);
	      gAlice->SetTrack(fTrackIt,-1,fIpart,p,origin,polar,0,kPPrimary,nt);
	  }
      }
  }
}









