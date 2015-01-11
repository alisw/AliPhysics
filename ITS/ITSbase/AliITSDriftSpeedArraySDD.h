#ifndef ALIITSDRIFTSPEEDARRAYSDD_H
#define ALIITSDRIFTSPEEDARRAYSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for a TOnjArray of the AliITSDriftSpeedSDD objects      //
// from 1 run (1 AliITSDriftSpeedSDD for  each injector trigger  //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include<TObject.h>
#include<TObjArray.h>

class AliITSDriftSpeedSDD;


class AliITSDriftSpeedArraySDD : public TObject{
 public:
  AliITSDriftSpeedArraySDD();
  AliITSDriftSpeedArraySDD(Int_t numEv);
  virtual ~AliITSDriftSpeedArraySDD() {};

  void AddDriftSpeed(AliITSDriftSpeedSDD* drSpeed);
  void SetInjectorStatus(UInt_t status=0x3E000000){fInjectorStatus=status;}
  void PrintAll() const;
  UInt_t GetTimestamp(Int_t iElement);
  UInt_t GetInjectorStatus() const {return fInjectorStatus;}
  Double_t GetDriftSpeed(Int_t iEvent, Double_t iAnode);
  AliITSDriftSpeedSDD* GetDriftSpeedObject(Int_t iEvent) const{
    if(iEvent>=0 && iEvent<fNEvents) return (AliITSDriftSpeedSDD*)fDriftSpeedSDD.At(iEvent);
    else return 0;
  }

 protected:  
  Int_t fNEvents;               // number of drift speed determination
  TObjArray fDriftSpeedSDD; // array of AliITSDriftSpeedSDD objects
  UInt_t fInjectorStatus;   // encoded info on injector status
                            // see AliITSOnlineSDDInjectors for definition

  ClassDef(AliITSDriftSpeedArraySDD,3);
};
#endif
