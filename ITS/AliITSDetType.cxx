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

////////////////////////////////////////////////////////////////////////
// This Class owns the classes needed to to detector simulations and
// reconstruction. This includes the detector segmentation classes,
// the detector responce classes, the detector simulatin classes, and
// the detector reconstruction (clustering) classes for all of the ITS
// detectors.
////////////////////////////////////////////////////////////////////////
#include "AliITSDetType.h"
#include "AliITSClusterFinder.h"
#include "AliITSsimulation.h"


ClassImp(AliITSDetType)	 

//______________________________________________________________________
AliITSDetType::AliITSDetType(){
    // constructor

    fSegmentation = 0;
    fResponse     = 0;
    fSimulation   = 0;
    fReconst      = 0;
}
//----------------------------------------------------------------------
AliITSDetType::~AliITSDetType(){
    // destructor

    if(fSegmentation!=0) delete fSegmentation; fSegmentation = 0;
    if(fResponse!=0)     delete fResponse;     fResponse     = 0;
    if(fSimulation!=0)   delete fSimulation;   fSimulation   = 0;
    if(fReconst!=0)      delete fReconst;      fReconst      = 0;
}
//______________________________________________________________________
AliITSDetType::AliITSDetType(const AliITSDetType &source) : TObject(source){
  //     Copy Constructor 

  if(&source == this) return;
  this->fReconst        = source.fReconst;
  this->fSimulation     = source.fSimulation;
  this->fResponse       = source.fResponse;
  this->fSegmentation   = source.fSegmentation;
  this->fDigClassName   = source.fDigClassName;
  this->fClustClassName = source.fClustClassName;
  return;
}
//______________________________________________________________________
AliITSDetType& AliITSDetType::operator=(const AliITSDetType &source){
  //    Assignment operator

  if(&source == this) return *this;
  this->fReconst        = source.fReconst;
  this->fSimulation     = source.fSimulation;
  this->fResponse       = source.fResponse;
  this->fSegmentation   = source.fSegmentation;
  this->fDigClassName   = source.fDigClassName;
  this->fClustClassName = source.fClustClassName;
  return *this;  
}
