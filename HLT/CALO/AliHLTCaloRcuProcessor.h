//-*- Mode: C++ -*-
// $Id: AliHLTCALORcuProcessor.h 29824 2008-11-10 13:43:55Z richterm $

#ifndef ALIHLTCALORCUPROCESSOR_H
#define ALIHLTCALORCUPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliHLTCaloProcessor.h"
//#include "AliHLTCALORcuProperties.h"


class AliHLTCaloRcuProcessor:public AliHLTCaloProcessor
{
public:
  AliHLTCaloRcuProcessor();
  virtual ~AliHLTCaloRcuProcessor();
  virtual AliHLTComponentDataType GetOutputDataType() =0;
  
private:
  AliHLTCaloRcuProcessor (const AliHLTCaloRcuProcessor & );
  AliHLTCaloRcuProcessor   & operator = (const  AliHLTCaloRcuProcessor  &);
    
};

#endif


