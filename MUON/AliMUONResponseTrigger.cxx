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
Revision 1.3  2000/07/03 11:54:57  morsch
AliMUONSegmentation and AliMUONHitMap have been replaced by AliSegmentation and AliHitMap in STEER
The methods GetPadIxy and GetPadXxy of AliMUONSegmentation have changed name to GetPadI and GetPadC.

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.1  2000/06/09 21:48:28  morsch
Code from AliMUONSegResTrigger.cxx

*/

#include "AliMUONResponseTrigger.h"
#include "AliSegmentation.h"
#include <TMath.h>
#include <TRandom.h>
#include <iostream.h> 

ClassImp(AliMUONResponseTrigger)

//------------------------------------------------------------------   
Int_t AliMUONResponseTrigger::SetGenerCluster(){
// nothing to be done except return 0
  return 0;
} 

//------------------------------------------------------------------   
Int_t AliMUONResponseTrigger::DigitResponse(Int_t digit){
//  only digital (0/1) information available
  if (digit) digit=1;
  return digit;
}






