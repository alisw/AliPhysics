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
  Revision 1.1  2000/06/12 15:28:06  jbarbosa
  Cleaned up version.

*/


#include "AliRICHRecHit.h"
 
ClassImp(AliRICHRecHit)

AliRICHRecHit::AliRICHRecHit(Int_t id, Float_t *rechit, Float_t *photons, Int_t *padsx, Int_t* padsy)
{
    //
    // Creates a RICH rec. hit object
    //
    fTheta        = rechit[0];
    fPhi          = rechit[1];
    fOmega        = rechit[2];
    fX            = rechit[3];
    fY            = rechit[4];
    fEmissPoint   = rechit[5];
    fGoodPhotons  = (Int_t)rechit[6];
    for(Int_t i=0; i<100; i++) {
	fCerPerPhoton[i]  = photons[i];
	fPadsUsedX[i]    = padsx[i];
	fPadsUsedY[i]    = padsy[i];
    }
    
}


