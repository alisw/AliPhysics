#ifndef ALIPHOSRAWDIGIPRODUCER_H
#define ALIPHOSRAWDIGIPRODUCER_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

//This class produces PHOS digits of one event
//using AliPHOSRawDecoder. See cxx source for use case.

class AliPHOSRawDecoder;

class AliPHOSRawDigiProducer {

public:

  AliPHOSRawDigiProducer() {}
  virtual ~AliPHOSRawDigiProducer() {}

  void MakeDigits(TClonesArray *digits, AliPHOSRawDecoder* decoder);

  ClassDef(AliPHOSRawDigiProducer,1)
};

#endif
