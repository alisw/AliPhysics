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

#include "AliMUONPadHit.h"

ClassImp(AliMUONPadHit)
 
//___________________________________________
AliMUONPadHit::AliMUONPadHit()
  :  TObject(),
     fHitNumber(0),
     fCathode(0),
     fQ(0),
     fPadX(0),
     fPadY(0),
     fQpad(0),
     fRSec(0)   
{
// Default constructor
}

//___________________________________________
AliMUONPadHit::AliMUONPadHit(Int_t *clhits)
  : TObject()
{
// Constructor
   fHitNumber=clhits[0];
   fCathode=clhits[1];
   fQ=clhits[2];
   fPadX=clhits[3];
   fPadY=clhits[4];
   fQpad=clhits[5];
   fRSec=clhits[6];
}

