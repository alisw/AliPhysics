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

// Read background particles from a boundary source
// Very specialized generator to simulate background from beam halo.
// The input file is a text file specially prepared 
// for this purpose.
// Author: andreas.morsch@cern.ch

#include <stdlib.h>

#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TSystem.h>

#include "AliGenHaloProtvino.h"
#include "AliRun.h"

ClassImp(AliGenHaloProtvino)

AliGenHaloProtvino::AliGenHaloProtvino()
    :AliGenerator(-1)
{
// Constructor
    
    fName  = "HaloProtvino";
    fTitle = "Halo from LHC Tunnel";
//
//  Read all particles
    fNpart = -1;
    fFile  =  0;
    fSide  =  1;
//
    SetRunPeriod();
    SetTimePerEvent();
    SetAnalog(0);
}

AliGenHaloProtvino::AliGenHaloProtvino(Int_t npart)
    :AliGenerator(npart)
{
// Constructor
    fName = "Halo";
    fTitle= "Halo from LHC Tunnel";
//
    fNpart   = npart;
    fFile    = 0;
    fSide    = 1;
//
    SetRunPeriod();
    SetTimePerEvent();
    SetAnalog(0);
}

AliGenHaloProtvino::AliGenHaloProtvino(const AliGenHaloProtvino & HaloProtvino):
    AliGenerator(HaloProtvino)
{
// Copy constructor
    HaloProtvino.Copy(*this);
}


//____________________________________________________________
AliGenHaloProtvino::~AliGenHaloProtvino()
{
// Destructor
}

//____________________________________________________________
void AliGenHaloProtvino::Init() 
{
// Initialisation
    fFile = fopen(fFileName,"r");
    if (fFile) {
	printf("\n File %s opened for reading, %p ! \n ",  fFileName.Data(), fFile);
    } else {
	printf("\n Opening of file %s failed,  %p ! \n ",  fFileName.Data(), fFile);
    }
//
//
//
//    Read file with gas pressure values
    char *name;
    
    name = gSystem->ExpandPathName("$(ALICE_ROOT)/LHC/gasPressure.dat" );
    FILE* file = fopen(name, "r");
    Float_t z;
    Int_t i, j;
    
//
//  Ring 1   
// 
    for (i = 0; i < 21; i++)
    {
	fscanf(file, "%f %f %f %f %f %f", &z, 
	       &fG1[i][0], &fG1[i][1], &fG1[i][2], &fG1[i][3], &fG1[i][4]);
	if (i > 0) {
	    fZ1[i] = fZ1[i-1] + z;
	} else {
	    fZ1[i] = 20.;
	}
  }
//
// Ring 2
//
    for (i = 0; i < 21; i++)
    {
	fscanf(file, "%f %f %f %f %f %f", &z, 
	       &fG2[i][0], &fG2[i][1], &fG2[i][2], &fG2[i][3], &fG2[i][4]);
	if (i > 0) {
	    fZ2[i] = fZ2[i-1] + z;
	} else {
	    fZ2[i] = 20.;
	}
    }
//
//  Transform into interaction rates
//
    const Float_t kCrossSection = 0.094e-28;     // m^2
    Float_t pFlux[5] = {0.2, 0.2, 0.3, 0.3, 1.0};

    for (j = 0; j <  5; j++) {
	pFlux[j] *= 1.e11/25.e-9;
	for (i = 0; i < 21; i++)  
	{
	    fG1[i][j] = fG1[i][j] * kCrossSection * pFlux[j]; // 1/m/s 
	    fG2[i][j] = fG2[i][j] * kCrossSection * pFlux[j]; // 1/m/s
	}
    }
    

    Float_t sum1 = 0.;
    Float_t sum2 = 0.;
    
    for (Int_t i = 0; i < 300; i++) {
	Float_t z = 20.+i*1.;
	z*=100;
	Float_t wgt1 = GassPressureWeight(z);
	Float_t wgt2 = GassPressureWeight(-z);
//	printf("weight: %f %f %f %f %f \n", z, wgt1, wgt2, fZ1[20], fZ2[20]);
	sum1 += wgt1;
	sum2 += wgt2;
    }
    sum1/=250.;
    sum2/=250.;
    printf("\n %f %f \n \n", sum1, sum2);
}

//____________________________________________________________
void AliGenHaloProtvino::Generate()
{
// Generate from input file
 
  Float_t polar[3]= {0,0,0};
  Float_t origin[3];
  Float_t p[3], p0;
  Float_t tz, txy;
  Float_t amass;
  //
  Int_t ncols, nt;
  Int_t nskip = 0;
  Int_t nread = 0;

  Float_t* zPrimary = new Float_t [fNpart];
  Int_t  * inuc     = new Int_t   [fNpart];
  Int_t  * ipart    = new Int_t   [fNpart];
  Float_t* wgt      = new Float_t [fNpart]; 
  Float_t* ekin     = new Float_t [fNpart];
  Float_t* vx       = new Float_t [fNpart];
  Float_t* vy       = new Float_t [fNpart];
  Float_t* tx       = new Float_t [fNpart];
  Float_t* ty       = new Float_t [fNpart];
  
  Float_t zVertexOld = -1.e10;
  Int_t   nInt       = 0;        // Counts number of interactions
  
  while(1) {
//
// Load event into array
//
      ncols = fscanf(fFile,"%f %d %d %f %f %f %f %f %f",
		     &zPrimary[nread], &inuc[nread], &ipart[nread], &wgt[nread], 
		     &ekin[nread], &vx[nread], &vy[nread],
		     &tx[nread], &ty[nread]);
      
      if (ncols < 0) break;
// Skip fNskip events
      nskip++;
      if (fNpart !=-1 && nskip <= fNskip) continue;
// Count interactions
      if (zPrimary[nread] != zVertexOld) {
	  nInt++;
	  zVertexOld = zPrimary[nread];
      }
// Count tracks      
      nread++;
      if (fNpart !=-1 && nread > fNpart) break;
  }
//
// Mean time between interactions
//
  Float_t dT = fTimePerEvent/nInt;   // sec 
  Float_t t  = 0;                    // sec
  
//
// Loop over primaries
//
  zVertexOld   = -1.e10;
  Double_t arg = 0.;
  
  for (Int_t nprim = 0; nprim < fNpart; nprim++) 
  {
      amass = TDatabasePDG::Instance()->GetParticle(ipart[nprim])->Mass();

      //
      // Momentum vector
      //
      p0=sqrt(ekin[nprim]*ekin[nprim] + 2.*amass*ekin[nprim]);
      
      txy=TMath::Sqrt(tx[nprim]*tx[nprim]+ty[nprim]*ty[nprim]);
      if (txy == 1.) {
	  tz=0;
      } else {
	  tz=-TMath::Sqrt(1.-txy);
      }
    
      p[0] = p0*tx[nprim];
      p[1] = p0*ty[nprim];
      p[2] =-p0*tz;
      
      origin[0] = vx[nprim];
      origin[1] = vy[nprim];
      origin[2] = -2196.5;

      //
      //
      // Particle weight

      Float_t originP[3] = {0., 0., 0.};
      originP[2] = zPrimary[nprim];
      
      Float_t pP[3] = {0., 0., 0.};
      Int_t ntP;
      
      if (fSide == -1) {
	  originP[2] = -zPrimary[nprim];
	  origin[2]  = -origin[2];
	  p[2]       = -p[2];
      }

      //
      // Time
      //
      if (zPrimary[nprim] != zVertexOld) {
	  while(arg==0.) arg = gRandom->Rndm();
	  t -= dT*TMath::Log(arg);              // (sec)   
	  zVertexOld = zPrimary[nprim];
      }

//    Get statistical weight according to local gas-pressure
//
      fParentWeight=wgt[nprim]*GassPressureWeight(zPrimary[nprim]);

      if (!fAnalog || gRandom->Rndm() < fParentWeight) {
//    Pass parent particle
//
	  PushTrack(0,-1,kProton,pP,originP,polar,t,kPNoProcess,ntP, fParentWeight);
	  KeepTrack(ntP);
	  PushTrack(fTrackIt,ntP,ipart[nprim],p,origin,polar,t,kPNoProcess,nt,fParentWeight);
      }

      //
      // Both sides are considered
      //

      if (fSide > 1) {
	  fParentWeight=wgt[nprim]*GassPressureWeight(-zPrimary[nprim]);
	  if (!fAnalog || gRandom->Rndm() < fParentWeight) {
	      origin[2]  = -origin[2];
	      originP[2] = -originP[2];
	      p[2]=-p[2];
	      PushTrack(0,-1,kProton,pP,originP,polar,t,kPNoProcess,ntP, fParentWeight);
	      KeepTrack(ntP);
	      PushTrack(fTrackIt,ntP,ipart[nprim],p,origin,polar,t,kPNoProcess,nt,fParentWeight);
	  }
      }
      SetHighWaterMark(nt);
  }
  delete [] zPrimary;
  delete [] inuc;    
  delete [] ipart;   
  delete [] wgt;     
  delete [] ekin;    
  delete [] vx;      
  delete [] vy;      
  delete [] tx;      
  delete [] ty;      
}
 

AliGenHaloProtvino& AliGenHaloProtvino::operator=(const  AliGenHaloProtvino& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}



Float_t AliGenHaloProtvino::GassPressureWeight(Float_t zPrimary)
{
//
// Return z-dependent gasspressure weight = interaction rate [1/m/s].
//
    Float_t weight = 0.;
    zPrimary/=100.;        // m
    Float_t zAbs = TMath::Abs(zPrimary);
    if (zPrimary > 0.) 
    {
	if (zAbs > fZ1[20]) {
	    weight = 2.e4;
	} else {
	    for (Int_t i = 1; i < 21; i++) {
		if (zAbs < fZ1[i]) {
		    weight = fG1[i][fRunPeriod];
		    break;
		}
	    }
	}
    } else {
	if (zAbs > fZ2[20]) {
	    weight = 2.e4;
	} else {
	    for (Int_t i = 1; i < 21; i++) {
		if (zAbs < fZ2[i]) {
		    weight = fG2[i][fRunPeriod];
		    break;
		}
	    }
	}
    }
    return weight;
}

void AliGenHaloProtvino::Copy(AliGenHaloProtvino&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}


/*
# Title:    README file for the sources of IR8 machine induced background
# Author:   Vadim Talanov <Vadim.Talanov@cern.ch>
# Modified: 12-12-2000 

0. Overview

	There are three files, named ring.one.beta.[01,10,50].m, which
	contain the lists of background particles, induced by proton losses
	upstream of IP8 in the LHC ring one, for the beta* values of 1, 10
	and 50 m, respectively.

1. File contents

	Each line in the files contains the coordinates of particle track
	crossing with the infinite plane, positioned at z=-1m, together with
	the physical properties of corresponding particle, namely:

	S  - S coordinate of the primary interaction vertex, cm;
	N  - type of the gas nuclei at interaction, 1 is H, 2 - C and 3 - O;
	I  - particle ID in PDG particle numbering scheme;
	W  - particle weight;
	E  - particle kinetic energy, GeV;
	X  - x coordinate of the crossing point, cm;
	Y  - y coordinate of the crossing point, cm;
	Dx - x direction cosine;
	Dy - y direction cosine.

2. Normalisation

	Each file is given per unity of linear density of proton inelastic
	interactions with the gas nuclei, [1 inelastic interaction/m].

# ~/vtalanov/public/README.mib: the end.

*/




