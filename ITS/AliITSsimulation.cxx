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

#include "AliITSsimulation.h"

ClassImp(AliITSsimulation)	

AliITSsimulation::AliITSsimulation() 
{
  // constructor
    fSegmentation=0;
    fResponse=0;
}

//__________________________________________________________________________
AliITSsimulation::AliITSsimulation(const AliITSsimulation &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fResponse = source.fResponse;
  this->fSegmentation = source.fSegmentation;
  return;
}

//_________________________________________________________________________
AliITSsimulation& 
  AliITSsimulation::operator=(const AliITSsimulation &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fResponse = source.fResponse; 
  this->fSegmentation = source.fSegmentation;
  return *this;
}
