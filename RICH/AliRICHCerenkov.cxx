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
  Revision 1.1  2000/06/12 15:16:59  jbarbosa
  Cleaned up version.

*/
#include "AliRICHCerenkov.h"

ClassImp(AliRICHCerenkov)
//___________________________________________
AliRICHCerenkov::AliRICHCerenkov(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
    AliHit(shunt, track)
{
// Constructor for object AliRICHCerenkov
    fChamber=vol[0];
    fX=hits[1];
    fY=hits[2];
    fZ=hits[3];
    fTheta=hits[4];
    fPhi=hits[5];
    fTlength=hits[6];
    fEloss=hits[7];
    fPHfirst=(Int_t) hits[8];
    fPHlast=(Int_t) hits[9];
    fCMother=Int_t(hits[10]);
    fIndex = hits[11];
    fProduction = hits[12];  
    fLoss=hits[13];
    fMomX=hits[14];
    fMomY=hits[15];
    fMomZ=hits[16];
    fNPads=hits[17];
    fCerenkovAngle=hits[18];
}

