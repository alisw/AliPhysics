//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************


/* $Id$ */

//
// BCM Hit Class
// 
// Detector Numberig 
// 11, 12, 13, 14 for z > 0
// 21, 22, 23, 24 for z < 0
//
// andreas.morsch@cern.ch


#include "AliBCMHit.h"
ClassImp(AliBCMHit)


AliBCMHit::AliBCMHit():
    AliHit(),
    fId(-1),
    fEdep(0.),
    fTime(0.)
{
    // Default constructor
}

AliBCMHit::AliBCMHit(Int_t shunt, Int_t track, Double_t x[4], Int_t isens, Float_t edep):
    AliHit(shunt, track),
    fId(isens),
    fEdep(edep),
    fTime(x[3])
{
    // Constructor
    fX = x[0];
    fY = x[1];
    fZ = x[2];
}
