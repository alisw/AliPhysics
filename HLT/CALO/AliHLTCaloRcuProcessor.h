//-*- Mode: C++ -*-
// $Id: AliHLTCaloRcuProcessor.h 29824 2008-11-10 13:43:55Z richterm $

#ifndef ALIHLTCALORCUPROCESSOR_H
#define ALIHLTCALORCUPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliHLTCaloProcessor.h"
//#include "AliHLTCaloRcuProperties.h"

//class  AliHLTCaloRcuProcessor : public AliHLTCaloProcessor, public AliHLTCaloRcuProperties
class  AliHLTCaloRcuProcessor : public AliHLTCaloProcessor
{
 public:
  AliHLTCaloRcuProcessor();
  virtual ~AliHLTCaloRcuProcessor();

 private:
  AliHLTCaloRcuProcessor (const AliHLTCaloRcuProcessor & );
  AliHLTCaloRcuProcessor   & operator = (const  AliHLTCaloRcuProcessor  &);



};

#endif


