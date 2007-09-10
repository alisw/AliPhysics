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

//-------------------------------------------------------------------------
//                      Implementation of   Class AliESDCaloTrigger
//   This is a class that summarizes the Trigger Data of EMCal and Phos
//   for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------


#include "AliESDCaloTrigger.h"

ClassImp(AliESDCaloTrigger)

AliESDCaloTrigger::AliESDCaloTrigger() : 
  TNamed(),
  fTriggerAmplitudes(0x0),
  fTriggerPosition(0x0)
{
}

AliESDCaloTrigger::AliESDCaloTrigger(const AliESDCaloTrigger &ctrig) : 
  TNamed(ctrig),
  fTriggerAmplitudes(ctrig.fTriggerAmplitudes),
  fTriggerPosition(ctrig.fTriggerPosition)
{
}

AliESDCaloTrigger::~AliESDCaloTrigger()
{
  delete fTriggerAmplitudes; fTriggerAmplitudes = 0;
  delete fTriggerPosition; fTriggerPosition = 0;
}

AliESDCaloTrigger& AliESDCaloTrigger::operator=(const AliESDCaloTrigger& ctrig)
{
  // assigment operator
  if(this!=&ctrig) {
    TNamed::operator=(ctrig);
    // CKB dont't want to create leak if fTriggerAmp points to 
    // something already 
    delete fTriggerAmplitudes;
    fTriggerAmplitudes = new TArrayF(*ctrig.fTriggerAmplitudes);
    delete fTriggerPosition;    
    fTriggerPosition = new TArrayF(*ctrig.fTriggerPosition);
  } 
  return *this;
}

void AliESDCaloTrigger::Reset()
{
  // simple reset
  if( fTriggerAmplitudes){  
    fTriggerAmplitudes->Reset();
  }
  if( fTriggerPosition){
    fTriggerPosition->Reset();
  }
}


