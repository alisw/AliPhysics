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

#include "AliMUONConstants.h"


ClassImp(AliMUONConstants)

Int_t   AliMUONConstants::fgNCh = 14;
Int_t   AliMUONConstants::fgNTrackingCh = 10;
Int_t   AliMUONConstants::fgNTriggerCh = 4;
Int_t   AliMUONConstants::fgNTriggerCircuit = 234;
Float_t AliMUONConstants::fgDefaultChamberZ[14] =
{518., 538., 680., 700., 965., 985., 1239., 1259., 1439., 1459.,
		   1610, 1625., 1710., 1725.}; 
Float_t  AliMUONConstants::fgDmin[7] = { 35.,  47.,  66.,   80.,  80., 100., 100.};    
Float_t  AliMUONConstants::fgDmax[7]  = {183., 245., 316.6,  520.,  520., 830., 880.};  

void AliMUONConstants::Streamer(TBuffer &R__b) {} 
