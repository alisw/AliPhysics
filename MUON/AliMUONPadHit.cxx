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
Revision 1.1.2.1  2000/06/09 22:02:45  morsch
Was before in DataStructures.cxx

*/

#include "AliMUONPadHit.h"

ClassImp(AliMUONPadHit)
 
//___________________________________________
AliMUONPadHit::AliMUONPadHit(Int_t *clhits)
{
   fHitNumber=clhits[0];
   fCathode=clhits[1];
   fQ=clhits[2];
   fPadX=clhits[3];
   fPadY=clhits[4];
   fQpad=clhits[5];
   fRSec=clhits[6];
}

