//-*- Mode: C++ -*-
// $Id: AliHLTCALORcuProcessor.h 29824 2008-11-10 13:43:55Z richterm $

#ifndef ALIHLTCALORCUPROCESSOR_H
#define ALIHLTCALORCUPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//#include "AliHLTCALOProcessor.h"
//#include "AliHLTCALORcuProperties.h"

//class  AliHLTCALORcuProcessor : public AliHLTCALOProcessor, public AliHLTCALORcuProperties
class  AliHLTCALORcuProcessor : public AliHLTCALOProcessor
{
 public:
  AliHLTCALORcuProcessor();
  virtual ~AliHLTCALORcuProcessor();

 private:
  AliHLTCALORcuProcessor (const AliHLTCALORcuProcessor & );
  AliHLTCALORcuProcessor   & operator = (const  AliHLTCALORcuProcessor  &);



};

#endif


