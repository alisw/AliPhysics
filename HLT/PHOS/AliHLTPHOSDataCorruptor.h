//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSDATACORRUPTOR_H
#define ALIHLTPHOSDATACORRUPTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include <Rtypes.h>

class AliHLTPHOSPulseGenerator;
class TRandom;

class AliHLTPHOSDataCorruptor
{
 public:
  AliHLTPHOSDataCorruptor();
  virtual ~ AliHLTPHOSDataCorruptor();
  AliHLTPHOSDataCorruptor(const AliHLTPHOSDataCorruptor & );
  AliHLTPHOSDataCorruptor & operator = (const AliHLTPHOSDataCorruptor &)
   {
      return *this;
   };

  AliHLTPHOSPulseGenerator *fPulseGeneratorPtr;
  TRandom *fRandomGeneratorPtr;

  void MakeCorruptedData(Double_t *dataArray, int N);
  void MakeCorruptedDataTest(Double_t *dataArray, int N);
  void FlipBit(int *sample, int n);
};



#endif
