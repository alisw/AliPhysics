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

#include "AliGenHalo.h"
#include "AliRun.h"

ClassImp(AliGenHalo)

AliGenHalo::AliGenHalo()
    :AliGenerator(-1),
     fp(0),
     fFileName(0)
{
// Constructor
    fName="Halo";
    fTitle="Halo from LHC Tunnel";
//
//  Read all particles
    fNpart=-1;
}

AliGenHalo::AliGenHalo(Int_t npart)
    :AliGenerator(npart),
     fp(0),
     fFileName(0)
{
// Constructor
    fName="Halo";
    fTitle="Halo from LHC Tunnel";
//
//  Read all particles
    fNpart=-1;
}

//____________________________________________________________
AliGenHalo::~AliGenHalo()
{
// Destructor
}

//____________________________________________________________
void AliGenHalo::Init() 
{
// Initialisation
}

//____________________________________________________________
void AliGenHalo::Generate()
{
// Generate from input file
    fp = fopen(fFileName,"r");
    if (fp) {
	printf("\n File %s opened for reading ! \n ", (char*) &fFileName);
    } else {
	printf("\n Opening of file %s failed ! \n ",  (char*) &fFileName);
    }
//
// MARS particle codes
  const Int_t kmars[12]={0,kProton,kNeutron,kPiPlus,kPiMinus,kKPlus,kKMinus,
			 kMuonPlus,kMuonMinus,kGamma,kElectron,kPositron};
 
  Float_t polar[3]= {0,0,0};
  Float_t origin[3];
  Float_t p[3], p0;
  Float_t ekin, wgt, tx, ty, tz, txy;
  Float_t amass;
  //
  Int_t ipart, ncols, nt;
  
  Int_t nread=0;
  origin[2]=2650;
  
  while(1) {
      ncols = fscanf(fp,"%i %f %f %f %f %f %f",
		     &ipart, &ekin, &wgt, 
		     &origin[0], &origin[1],
		     &tx, &ty);
      if (ncols < 0) break;
      nread++;
      if (fNpart !=-1 && nread > fNpart) break;
      ipart = kmars[ipart];
      amass = TDatabasePDG::Instance()->GetParticle(ipart)->Mass();
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
      fParentWeight=wgt;
      PushTrack(fTrackIt,-1,ipart,p,origin,polar,0,kPNoProcess,nt,fParentWeight);
//    PushTrack(fTrackIt,-1,ipart,p,origin,polar,0,"Halo+",nt,fParentWeight);
      origin[2]=-origin[2];
      p[2]=-p[2];
      PushTrack(fTrackIt,-1,ipart,p,origin,polar,0,kPNoProcess,nt,fParentWeight);
//    PushTrack(fTrackIt,-1,ipart,p,origin,polar,0,"Halo-",nt,fParentWeight);
      origin[2]=-origin[2];
      p[2]=-p[2];
  }
}
 







