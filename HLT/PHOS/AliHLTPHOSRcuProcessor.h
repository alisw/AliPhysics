//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRCUPROCESSOR_H
#define ALIHLTPHOSRCUPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliHLTCaloProcessor.h"
//#include "AliHLTPHOSRcuProperties.h"

//class  AliHLTPHOSRcuProcessor : public AliHLTCaloProcessor, public AliHLTPHOSRcuProperties
class  AliHLTPHOSRcuProcessor : public AliHLTCaloProcessor

{
 public:
  AliHLTPHOSRcuProcessor();
  virtual ~AliHLTPHOSRcuProcessor();

 private:
  AliHLTPHOSRcuProcessor (const AliHLTPHOSRcuProcessor & );
  AliHLTPHOSRcuProcessor   & operator = (const  AliHLTPHOSRcuProcessor  &);



};

#endif


