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
Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.1  2000/06/09 22:02:14  morsch
Was before in DataStructures.cxx

*/
#include "AliMUONHit.h"

ClassImp(AliMUONHit)
 
//___________________________________________
    AliMUONHit::AliMUONHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
	AliHit(shunt, track)
{
// Constructor
    fChamber   = vol[0];
    fParticle  = hits[0];
    fX         = hits[1];
    fY         = hits[2];
    fZ         = hits[3];
    fTheta     = hits[4];
    fPhi       = hits[5];
    fTlength   = hits[6];
    fEloss     = hits[7];
    fPHfirst   = (Int_t) hits[8];
    fPHlast    = (Int_t) hits[9];
    fPTot      = hits[10];
    fPx        = hits[11];
    fPy        = hits[12];
    fPz        = hits[13];
    fAge       = hits[14];
}







