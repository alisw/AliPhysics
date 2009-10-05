//-*- Mode: C++ -*-
// $Id: AliHLTPHOSRawAnalyzerCrudeComponent.h 29824 2008-11-10 13:43:55Z richterm $

#ifndef ALIHLTCALORAWANALYZERCRUDECOMPONENT_H
#define ALIHLTCALORAWANALYZERCRUDECOMPONENT_H

#include "AliHLTCaloRawAnalyzerComponentv3.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 



class AliHLTCaloRawAnalyzerCrudeComponent: public AliHLTCaloRawAnalyzerComponentv3
{
 public:
  AliHLTCaloRawAnalyzerCrudeComponent();
  virtual ~AliHLTCaloRawAnalyzerCrudeComponent();
  AliHLTCaloRawAnalyzerCrudeComponent(const AliHLTCaloRawAnalyzerCrudeComponent & );
  // AliHLTCaloRawAnalyzerCrudeComponent & operator = (const AliHLTCaloRawAnalyzerCrudeComponent)
  //  {
  //   return *this;
  // };

  virtual int Deinit();
  virtual const char* GetComponentID();
  virtual AliHLTComponent* Spawn() = 0;
};

#endif
