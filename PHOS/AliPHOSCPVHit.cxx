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
  Revision 1.1  2000/11/03 16:49:53  schutz
  New class AliPHOSCPVHit

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

AliPHOSCPVHit::AliPHOSCPVHit(Int_t shunt, Int_t track, TLorentzVector p, Float_t *xy, Int_t ipart)
             : AliHit(shunt, track)
{
  //
  // Create a CPV hit object
  //

  fMomentum  = p;
  fX         = xy[0];
  fY         = xy[1];
  fIpart     = ipart;
}

//______________________________________________________________________________
void AliPHOSCPVHit::Print()
{
  //
  // Print CPV hit
  //

  printf("CPV hit: p  = (% .4f, % .4f, % .4f, % .4f) GeV,\n",
	fMomentum.Px(),fMomentum.Py(),fMomentum.Pz(),fMomentum.E());
  printf("         xy = (%8.4f, %8.4f) cm, ipart = %d\n",
	 fX,fY,fIpart);
}
