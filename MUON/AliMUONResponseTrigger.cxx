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


#include "AliMUONResponseTrigger.h"

ClassImp(AliMUONResponseTrigger)

//------------------------------------------------------------------   
AliMUONResponseTrigger::AliMUONResponseTrigger()
  : AliMUONResponse()
{
// Default constructor
}

//------------------------------------------------------------------   
Int_t AliMUONResponseTrigger::SetGenerCluster(){
// nothing to be done except return 0
  return 0; 
} 

//------------------------------------------------------------------   
Int_t AliMUONResponseTrigger::DigitResponse(Int_t digit, 
					    AliMUONTransientDigit* /*where*/){
//  only digital (0/1) information available
  if (digit) digit=1;
  return digit;
}






