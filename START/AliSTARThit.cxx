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
Revision 1.1.1.1  2000/02/15 15:20:19  kir
Initial import to mni CVS

*/ 
#include "AliSTARThit.h"

ClassImp(AliSTARThit)

AliSTARThit::AliSTARThit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
//Normal START hit ctor
  
   fVolume = vol[0];
   fPmt=vol[1];
   fX=hits[0];
   fY=hits[1];
   fZ=hits[2];
   fEdep=hits[3];
   fEtot=hits[4];
   fParticle=Int_t (hits[5]);
   fTime=hits[6];
}
 
