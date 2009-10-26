//-*- Mode: C++ -*-
// $Id: AliHLTCaloRawAnalyzerCrude.h 34264 2009-08-14 18:29:23Z odjuvsla $

#ifndef ALIHLTCALORAWANALYZERCRUDE_H
#define ALIHLTCALORAWANALYZERCRUDE_H
#include "AliHLTCaloRawAnalyzer.h"


/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


class AliHLTCaloRawAnalyzerCrude : public AliHLTCaloRawAnalyzer
{
 public:
  AliHLTCaloRawAnalyzerCrude();
  virtual ~AliHLTCaloRawAnalyzerCrude();
  virtual void Evaluate(int start = 0, int lenght = 100);
private:
};

#endif
