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

////////////////////////////////////////////////
//  Hit class CPV                             //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 28 September 2000          //
////////////////////////////////////////////////
 
// --- ROOT system ---

// --- Standard library ---
#include <stdio.h>

// --- galice header files ---
#include "AliPHOSCPVHit.h"


ClassImp(AliPHOSCPVHit)

//______________________________________________________________________________

AliPHOSCPVHit::AliPHOSCPVHit(TLorentzVector p, Float_t *xy, Int_t ipart)
{
  //
  // Create a CPV hit object
  //

  fMomentum  = p;
  fXhit      = xy[0];
  fYhit      = xy[1];
  fIpart     = ipart;
}

//______________________________________________________________________________
void AliPHOSCPVHit::Print()
{
  //
  // Print CPV hit
  //

  printf("CPV hit: p  = (% .4f, % .4f, % .4f, % .4f) GeV,\n",
	GetMomentum().Px(),GetMomentum().Py(),GetMomentum().Pz(),GetMomentum().E());
  printf("         xy = (%8.4f, %8.4f) cm, ipart = %d\n",
	 fXhit,fYhit,fIpart);
}
