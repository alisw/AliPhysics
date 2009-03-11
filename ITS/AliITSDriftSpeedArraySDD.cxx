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

#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSDriftSpeedSDD.h"
#include "AliLog.h"

ClassImp(AliITSDriftSpeedArraySDD)
//______________________________________________________________________
AliITSDriftSpeedArraySDD::AliITSDriftSpeedArraySDD():
TObject(),
fNEvents(0),
fDriftSpeedSDD(10){
  // default constructor
}
//______________________________________________________________________
AliITSDriftSpeedArraySDD::AliITSDriftSpeedArraySDD(Int_t numEv):
TObject(),
fNEvents(0),
fDriftSpeedSDD(numEv){
  // standard constructor
}
//______________________________________________________________________
void AliITSDriftSpeedArraySDD::AddDriftSpeed(AliITSDriftSpeedSDD* drSpeed){
  // adds an AliITSDriftSpeedSDD object in the array
  fDriftSpeedSDD.AddLast(drSpeed);
  fNEvents++;
}
//______________________________________________________________________
void AliITSDriftSpeedArraySDD::PrintAll() const{
  // print drift speed parameters for all elements in the array
  printf("Array Size=%d\n",fDriftSpeedSDD.GetSize());
  printf("Array Elements =%d\n",fNEvents);
  for(Int_t i=0;i<fNEvents; i++){
    printf("     ====== Array el. #%d ======\n",i);
    AliITSDriftSpeedSDD *d=(AliITSDriftSpeedSDD*)fDriftSpeedSDD.At(i);
    if(d) d->PrintDriftSpeedParameters();
  }
}
//______________________________________________________________________
Double_t AliITSDriftSpeedArraySDD::GetDriftSpeed(Int_t iEvent, Double_t iAnode){
  // returns drift speed for given event number and anode
  if(!fDriftSpeedSDD.IsSorted()) fDriftSpeedSDD.Sort();
  if(fNEvents==1){
    AliITSDriftSpeedSDD *d=(AliITSDriftSpeedSDD*)fDriftSpeedSDD.At(0);
    return d->GetDriftSpeedAtAnode(iAnode);
  }else{
    Int_t nInjEv1=-1;
    Int_t nInjEv2=-1;
    AliITSDriftSpeedSDD *d1=NULL;
    AliITSDriftSpeedSDD *d2=NULL;
    for(Int_t i=0;i<fNEvents; i++){
      d1=d2;
      d2=(AliITSDriftSpeedSDD*)fDriftSpeedSDD.At(i);
      nInjEv1=nInjEv2;
      if(d2!=0){
	nInjEv2=d2->GetEventNumber();
	if(nInjEv2>=iEvent){
	  if(i==0) d1=(AliITSDriftSpeedSDD*)fDriftSpeedSDD.At(i+1);
	  nInjEv1=d1->GetEventNumber();
	  break;
	}
      }
    }
    if(nInjEv1>=0 && nInjEv2>=0 && nInjEv1!=nInjEv2){
      Double_t v1=d1->GetDriftSpeedAtAnode(iAnode);
      Double_t v2=d2->GetDriftSpeedAtAnode(iAnode);
      Double_t vdrift=(v2-v1)*(iEvent-nInjEv1)/(nInjEv2-nInjEv1)+v1;
      return vdrift;
    }
  }
  AliWarning("Vdrift interpolation error\n");
  return -999.;
}
