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
 Revision 1.5  2001/05/31 06:58:38  fca
 Patch problem with destructor

 Revision 1.4  2001/05/01 14:47:45  nilsen
 Fixed destructor so that it destroyes the pointers fSegmentation, fResponse,
 fSimulation, and fReconst if they have been allocated. The two TStrings
 fDigClassName and fClustClassName shoud be destroyed automaticaly. This should
 fix a small memory leak associated with digitization and reconstruction.

*/

#include "AliITSDetType.h"
#include "AliITSClusterFinder.h"
#include "AliITSsimulation.h"

ClassImp(AliITSDetType)	 

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
  /* PH Commented because AliITSDetType is not the owner
    if(fSegmentation!=0) delete fSegmentation; fSegmentation = 0;
    //    if(fResponse!=0)     delete fResponse;     fResponse     = 0;
    if(fSimulation!=0)   delete fSimulation;   fSimulation   = 0;
    if(fReconst!=0)      delete fReconst;      fReconst      = 0;
  */
}
//______________________________________________________________________
AliITSDetType::AliITSDetType(const AliITSDetType &source){
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
