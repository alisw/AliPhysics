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

//-----------------------------------------------------------------
//
//   Virtual class to access calorimeter 
//   (EMCAL, PHOS, PMD, FMD) trigger data
//   Author: Salvatore Aiola
//
//-----------------------------------------------------------------

#include "AliVCaloTrigger.h"

ClassImp(AliVCaloTrigger)
  
  AliVCaloTrigger::AliVCaloTrigger(const AliVCaloTrigger& vtrg) :
    TNamed(vtrg) { ; } // Copy constructor

AliVCaloTrigger& AliVCaloTrigger::operator=(const AliVCaloTrigger& vtrg)
{ 
  //Assignment operator
  if (this!=&vtrg) { 
    TObject::operator=(vtrg); 
  }
  
  return *this; 
}


void AliVCaloTrigger::Copy(TObject &obj) const 
{	
  TNamed::Copy(obj);
}
