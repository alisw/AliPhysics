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

#include "AliITSDetType.h"

#include "AliITSClusterFinder.h"
#include "AliITSsimulation.h"

ClassImp(AliITSDetType)	 

AliITSDetType::AliITSDetType() 
{
  // constructor
    fSegmentation=0;
    fResponse=0;
    fSimulation=0;
    fReconst=0;
}

AliITSDetType::~AliITSDetType() 
{
  // destructor
}
//__________________________________________________________________________
AliITSDetType::AliITSDetType(const AliITSDetType &source){
  //     Copy Constructor 

  if(&source == this) return;
  this->fReconst = source.fReconst;
  this->fSimulation = source.fSimulation;
  this->fResponse = source.fResponse;
  this->fSegmentation = source.fSegmentation;
  this->fDigClassName = source.fDigClassName;
  this->fClustClassName = source.fClustClassName;
  return;
}

//_________________________________________________________________________
AliITSDetType& AliITSDetType::operator=(const AliITSDetType &source){
  //    Assignment operator
  if(&source == this) return *this;
  this->fReconst = source.fReconst;
  this->fSimulation = source.fSimulation;
  this->fResponse = source.fResponse;
  this->fSegmentation = source.fSegmentation;
  this->fDigClassName = source.fDigClassName;
  this->fClustClassName = source.fClustClassName;
  return *this;
  
}
