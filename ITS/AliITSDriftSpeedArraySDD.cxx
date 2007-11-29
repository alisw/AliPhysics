/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class with array of AliITSDriftSpeedSDD //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSDriftSpeedSDD.h"
#include "AliLog.h"

ClassImp(AliITSDriftSpeedArraySDD)
//______________________________________________________________________
AliITSDriftSpeedArraySDD::AliITSDriftSpeedArraySDD():
TObject(),
fNEvents(0),
fDriftSpeedSDD(0){
  // default constructor
  fDriftSpeedSDD=new TClonesArray("AliITSDriftSpeedSDD",100);
}
//______________________________________________________________________
AliITSDriftSpeedArraySDD::AliITSDriftSpeedArraySDD(Int_t numEv):
TObject(),
fNEvents(0),
fDriftSpeedSDD(0){
  // standard constructor
  fDriftSpeedSDD=new TClonesArray("AliITSDriftSpeedSDD",numEv);
}
//______________________________________________________________________
AliITSDriftSpeedArraySDD::AliITSDriftSpeedArraySDD(const AliITSDriftSpeedArraySDD& array):
TObject(),
fNEvents(array.fNEvents),
fDriftSpeedSDD(array.fDriftSpeedSDD){
  // copy constructor
}
//______________________________________________________________________
AliITSDriftSpeedArraySDD& AliITSDriftSpeedArraySDD::operator=(const AliITSDriftSpeedArraySDD& array){
  // assignment operator
  this->~AliITSDriftSpeedArraySDD();
  new(this) AliITSDriftSpeedArraySDD(array);
  return *this;
}

//______________________________________________________________________
AliITSDriftSpeedArraySDD::~AliITSDriftSpeedArraySDD(){
  // destructor
  if(fDriftSpeedSDD){
    fDriftSpeedSDD->Delete();
    delete fDriftSpeedSDD;
  }
}
//______________________________________________________________________
void AliITSDriftSpeedArraySDD::AddDriftSpeed(AliITSDriftSpeedSDD* drSpeed){
  // adds an AliITSDriftSpeedSDD object in the array
  TClonesArray &arr = *fDriftSpeedSDD;
  new(arr[fNEvents]) AliITSDriftSpeedSDD(*drSpeed);
  fNEvents++;
  
}
//______________________________________________________________________
void AliITSDriftSpeedArraySDD::PrintAll() const{
  // print drift speed parameters for all elements in the array
  printf("Array Size=%d\n",fDriftSpeedSDD->GetSize());
  printf("Array Elements =%d\n",fNEvents);
  for(Int_t i=0;i<fNEvents; i++){
    printf("     ====== Array el. #%d ======\n",i);
    AliITSDriftSpeedSDD *d=(AliITSDriftSpeedSDD*)fDriftSpeedSDD->At(i);
    if(d) d->PrintDriftSpeedParameters();
  }
}
//______________________________________________________________________
Float_t AliITSDriftSpeedArraySDD::GetDriftSpeed(Int_t iEvent, Float_t iAnode) const{
  // returns drift speed for given event number and anode
  if(!fDriftSpeedSDD->IsSorted()) fDriftSpeedSDD->Sort();
  if(fNEvents==1){
    AliITSDriftSpeedSDD *d=(AliITSDriftSpeedSDD*)fDriftSpeedSDD->At(0);
    return d->GetDriftSpeedAtAnode(iAnode);
  }else{
    Int_t nInjEv1=-1;
    Int_t nInjEv2=-1;
    AliITSDriftSpeedSDD *d1=NULL;
    AliITSDriftSpeedSDD *d2=NULL;
    for(Int_t i=0;i<fNEvents; i++){
      d1=d2;
      d2=(AliITSDriftSpeedSDD*)fDriftSpeedSDD->At(i);
      nInjEv1=nInjEv2;
      if(d2!=0){
	nInjEv2=d2->GetEventNumber();
	if(nInjEv2>=iEvent){
	  if(i==0) d1=(AliITSDriftSpeedSDD*)fDriftSpeedSDD->At(i+1);
	  nInjEv1=d1->GetEventNumber();
	  break;
	}
      }
    }
    if(nInjEv1>=0 && nInjEv2>=0 && nInjEv1!=nInjEv2){
      Float_t v1=d1->GetDriftSpeedAtAnode(iAnode);
      Float_t v2=d2->GetDriftSpeedAtAnode(iAnode);
      Float_t vdrift=(v2-v1)*(iEvent-nInjEv1)/(nInjEv2-nInjEv1)+v1;
      return vdrift;
    }
  }
  AliWarning("Vdrift interpolation error\n");
  return -999.;
}
