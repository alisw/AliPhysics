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

/* $Id$ */

// Realisation of AliGenerator that generates particles with
// vertices on a user defined grid.
// The vertex positions can be smeared. 
// Momentum vectors are defined through the methods provided by AliGenerator.
// Author: andreas.morsch@cern.ch

#include "AliGenScan.h"

 ClassImp(AliGenScan)
    
 AliGenScan::AliGenScan()
     :AliGenerator(-1), 
      fXCmin(0),
      fXCmax(0),
      fNx(1),
      fYCmin(0),
      fYCmax(0),
      fNy(1),
      fZmin(0),
      fZmax(0),
      fNz(1),
      fIpart(0)
{
// Constructor
//
//  Read all particles
    fNpart=-1;
}

AliGenScan::AliGenScan(Int_t npart)
    :AliGenerator(npart), 
      fXCmin(0),
      fXCmax(0),
      fNx(1),
      fYCmin(0),
      fYCmax(0),
      fNy(1),
      fZmin(0),
      fZmax(0),
      fNz(1),
      fIpart(0)
{
// Constructor
    fName  = "Scan";
    fTitle = "Generator for particles on a grid";
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
    fXCmin=xmin;
    fXCmax=xmax;
    fNx=nx;
    fYCmin=ymin;
    fYCmax=ymax;
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
  if (fNx > 0) {
      dx=(fXCmax-fXCmin)/fNx;
  } else {
      dx=1e10;
  }

  if (fNy > 0) {
      dy=(fYCmax-fYCmin)/fNy;
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
	      origin[0]=fXCmin+ix*dx+2*(random[0]-0.5)*fOsigma[0];
	      origin[1]=fYCmin+iy*dy+2*(random[1]-0.5)*fOsigma[1];
	      origin[2]=fZmin+iz*dz+2*(random[2]-0.5)*fOsigma[2];	     
	      pmom=fPMin+random[3]*(fPMax-fPMin);
	      theta=fThetaMin+random[4]*(fThetaMax-fThetaMin);
	      phi=fPhiMin+random[5]*(fPhiMax-fPhiMin);
	      p[0] = pmom*TMath::Cos(phi)*TMath::Sin(theta);
	      p[1] = pmom*TMath::Sin(phi)*TMath::Sin(theta);
	      p[2] = pmom*TMath::Cos(theta);
	      PushTrack(fTrackIt,-1,fIpart,p,origin,polar,0,kPPrimary,nt);
	  }
      }
  }
}









