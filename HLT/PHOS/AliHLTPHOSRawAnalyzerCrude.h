//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRAWANALYZERCRUDE_H
#define ALIHLTPHOSRAWANALYZERCRUDE_H
#include "AliHLTPHOSRawAnalyzer.h"


/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


class AliHLTPHOSRawAnalyzerCrude : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerCrude();
  virtual ~AliHLTPHOSRawAnalyzerCrude();
  virtual void Evaluate(int start = 0, int lenght = 100);
  // private:
  //  ClassDef(AliHLTPHOSRawAnalyzerCrude, 2) 
  
};

#endif
