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
*/

// Read background particles from a boundary source
// Very specialized generator to simulate background from beam halo.
// The input file is a text file specially prepared 
// for this purpose.
// Author: andreas.morsch@cern.ch

#include "AliGenHaloProtvino.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliPDG.h"

#include <TDatabasePDG.h>
#include <stdlib.h>

 ClassImp(AliGenHaloProtvino)
     AliGenHaloProtvino::AliGenHaloProtvino()
	 :AliGenerator(-1)
{
// Constructor
    fName="HaloProtvino";
    fTitle="Halo from LHC Tunnel";
//
//  Read all particles
    fNpart=-1;
    fp=0;
}

AliGenHaloProtvino::AliGenHaloProtvino(Int_t npart)
    :AliGenerator(npart)
{
// Constructor
    fName="Halo";
    fTitle="Halo from LHC Tunnel";
//
//  Read all particles
    fNpart=-1;
    fp=0;
}

AliGenHaloProtvino::AliGenHaloProtvino(const AliGenHaloProtvino & HaloProtvino)
{
// copy constructor
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
}

//____________________________________________________________
void AliGenHaloProtvino::Generate()
{
// Generate from input file
    FILE *fp = fopen(fFileName,"r");
    if (fp) {
	printf("\n File %s opened for reading ! \n ", (char*) &fFileName);
    } else {
	printf("\n Opening of file %s failed ! \n ",  (char*) &fFileName);
    }
 
  Float_t polar[3]= {0,0,0};
  Float_t origin[3];
  Float_t p[3], p0;
  Float_t ekin, wgt, tx, ty, tz, txy;
  Float_t zPrimary;
  Float_t amass;
  Int_t   inuc;
  //
  Int_t ipart, ncols, nt;
  
  Int_t nread=0;
  origin[2]=2650;
  
  while(1) {
    
      ncols = fscanf(fp,"%f %d %d %f %f %f %f %f %f",
		     &zPrimary, &inuc, &ipart, &wgt, 
		     &ekin, &origin[0], &origin[1],
		     &tx, &ty);
      printf(" \n %f %d %d %f %f %f %f %f %f",
		     zPrimary, inuc, ipart, wgt, 
		     ekin, origin[0], origin[1],
		     tx, ty);
 
      if (ncols < 0) break;
      nread++;
      if (fNpart !=-1 && nread > fNpart) break;



      amass = TDatabasePDG::Instance()->GetParticle(ipart)->Mass();

      //
      // Momentum vector
      //
      p0=sqrt(ekin*ekin + 2.*amass);
      
      txy=TMath::Sqrt(tx*tx+ty*ty);
      if (txy == 1.) {
	  tz=0;
      } else {
	  tz=-TMath::Sqrt(1.-txy);
      }
    
      p[0]=p0*tx;
      p[1]=p0*ty;
      p[2]=p0*tz;
      //
      // Origin: backtracking to tunnel entrance 
      // 
      Float_t zTunnelEntrance = -20.;
      origin[2] = zTunnelEntrance;
      Float_t dzBack = -1-zTunnelEntrance;
      Float_t dsBack = -dzBack/tz;
      
      origin[0]+=tx*dsBack;
      origin[1]+=ty*dsBack;

      //
      //
      // Particle weight


      fParentWeight=wgt*GassPressureWeight(zPrimary);
      gAlice->SetTrack(fTrackIt,-1,ipart,p,origin,polar,0,kPNoProcess,nt,fParentWeight);
      //
      // Assume particles come from two directions with same probability
      origin[2]=-origin[2];
      p[2]=-p[2];
      fParentWeight=wgt*GassPressureWeight(-zPrimary);
      gAlice->SetTrack(fTrackIt,-1,ipart,p,origin,polar,0,kPNoProcess,nt,fParentWeight);
      origin[2]=-origin[2];
      p[2]=-p[2];
  }
}
 

AliGenHaloProtvino& AliGenHaloProtvino::operator=(const  AliGenHaloProtvino& rhs)
{
// Assignment operator
    return *this;
}


Float_t AliGenHaloProtvino::GassPressureWeight(Float_t zPrimary)
{
  // Return z-dependent gasspressure weight
  //
  return 1.;
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




