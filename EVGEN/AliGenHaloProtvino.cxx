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
Revision 1.3  2001/07/27 17:09:36  morsch
Use local SetTrack, KeepTrack and SetHighWaterMark methods
to delegate either to local stack or to stack owned by AliRun.
(Piotr Skowronski, A.M.)

Revision 1.2  2001/06/14 12:15:27  morsch
Bugs corrected. SetSide() method added.

Revision 1.1  2001/01/23 15:04:33  morsch
Generator to read beam halo file from Protvino group.

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
    printf("\n Calling Default Constructor");
    
    fName  = "HaloProtvino";
    fTitle = "Halo from LHC Tunnel";
//
//  Read all particles
    fNpart = -1;
    fFile  =  0;
    fSide  =  1;
}

AliGenHaloProtvino::AliGenHaloProtvino(Int_t npart)
    :AliGenerator(npart)
{
// Constructor
    printf("\n Calling Constructor");
    fName = "Halo";
    fTitle= "Halo from LHC Tunnel";
//
    fNpart   = npart;
    fFile    = 0;
    fSide    = 1;
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
    fFile = fopen(fFileName,"r");
    if (fFile) {
	printf("\n File %s opened for reading, %p ! \n ",  fFileName.Data(), fFile);
    } else {
	printf("\n Opening of file %s failed,  %p ! \n ",  fFileName.Data(), fFile);
    }
}

//____________________________________________________________
void AliGenHaloProtvino::Generate()
{
// Generate from input file
 
  Float_t polar[3]= {0,0,0};
  Float_t origin[3];
  Float_t p[3], p0;
  Float_t ekin, wgt, tx, ty, tz, txy;
  Float_t zPrimary;
  Float_t amass;
  Int_t   inuc;
  //
  Int_t ipart, ncols, nt;
  Int_t nskip = 0;
  Int_t nread = 0;
  while(1) {
      ncols = fscanf(fFile,"%f %d %d %f %f %f %f %f %f",
		     &zPrimary, &inuc, &ipart, &wgt, 
		     &ekin, &origin[0], &origin[1],
		     &tx, &ty);
/*
      printf(" \n %f %d %d %f %f %f %f %f %f",
             zPrimary, inuc, ipart, wgt, 
             ekin, origin[0], origin[1],
             tx, ty);
*/
      
      if (ncols < 0) break;

      nskip++;
      if (fNpart !=-1 && nskip <= fNskip) continue;
      
      nread++;
      if (fNpart !=-1 && nread > fNpart) break;

      amass = TDatabasePDG::Instance()->GetParticle(ipart)->Mass();

      //
      // Momentum vector
      //
      p0=sqrt(ekin*ekin + 2.*amass*ekin);
      
      txy=TMath::Sqrt(tx*tx+ty*ty);
      if (txy == 1.) {
	  tz=0;
      } else {
	  tz=-TMath::Sqrt(1.-txy);
      }
    
      p[0]=p0*tx;
      p[1]=p0*ty;
      p[2]=-p0*tz;
      
      origin[2] = -2196.5;

      //
      //
      // Particle weight

      Float_t originP[3] = {0., 0., 0.};
      originP[2] = zPrimary;
      
      Float_t pP[3] = {0., 0., 0.};
      Int_t ntP;
      
      if (fSide == -1) {
	  originP[2] = -zPrimary;
	  origin[2]  = -origin[2];
	  p[2]       = -p[2];
      }

      SetTrack(0,-1,kProton,pP,originP,polar,0,kPNoProcess,ntP);
      KeepTrack(ntP);
      fParentWeight=wgt*GassPressureWeight(zPrimary);
      SetTrack(fTrackIt,ntP,ipart,p,origin,polar,0,kPNoProcess,nt,fParentWeight);
      SetHighWaterMark(nt);
      
      //
      // Assume particles come from two directions with same probability

      // Both Side are considered
      if (fSide > 1) {
          origin[2]=-origin[2];
          p[2]=-p[2];
          fParentWeight=wgt*GassPressureWeight(-zPrimary);
          SetTrack(fTrackIt,ntP,ipart,p,origin,polar,0,kPNoProcess,nt,fParentWeight);
      }
      
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
    Float_t weight = 500.;
    
    if (zPrimary > 45000.) weight = 2.e4;
    
  return weight;
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




