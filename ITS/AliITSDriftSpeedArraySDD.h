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

  void PrintAll() const;
  Double_t GetDriftSpeed(Int_t iEvent, Double_t iAnode);

 protected:  
  Int_t fNEvents;               // number of drift speed determination
  TObjArray fDriftSpeedSDD; // array of AliITSDriftSpeedSDD objects
  ClassDef(AliITSDriftSpeedArraySDD,2);
};
#endif
