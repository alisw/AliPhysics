#include "AliGenScan.h"
#include <stdlib.h>
#include "AliRun.h"
 ClassImp(AliGenScan)
    
 AliGenScan::AliGenScan()
	 :AliGenerator(-1)
{
//
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
{}

void AliGenScan::SetRange(Int_t nx, Float_t xmin, Float_t xmax,
		     Int_t ny, Float_t ymin, Float_t ymax,
		     Int_t nz, Float_t zmin, Float_t zmax)
{
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
  AliMC* pMC = AliMC::GetMC();
  
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
	      pMC->Rndm(random,6);
	      origin[0]=fXmin+ix*dx+2*(random[0]-0.5)*fOsigma[0];
	      origin[1]=fYmin+iy*dy+2*(random[1]-0.5)*fOsigma[1];
	      origin[2]=fZmin+iz*dz+2*(random[2]-0.5)*fOsigma[2];	     
	      pmom=fPMin+random[3]*(fPMax-fPMin);
	      theta=fThetaMin+random[4]*(fThetaMax-fThetaMin);
	      phi=fPhiMin+random[5]*(fPhiMax-fPhiMin);
	      p[0] = pmom*TMath::Cos(phi)*TMath::Sin(theta);
	      p[1] = pmom*TMath::Sin(phi)*TMath::Sin(theta);
	      p[2] = pmom*TMath::Cos(theta);
	      gAlice->SetTrack(1,-1,fIpart,p,origin,polar,0,"Primary",nt);
	  }
      }
  }
}









