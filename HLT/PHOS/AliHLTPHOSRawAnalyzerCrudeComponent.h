//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRAWANALYZERCRUDECOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCRUDECOMPONENT_H

//#include "AliHLTPHOSRawAnalyzerComponent.h"
#include "AliHLTPHOSRawAnalyzerComponentv3.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 



//class AliHLTPHOSRawAnalyzerCrudeComponent: public AliHLTPHOSRawAnalyzerComponent
class AliHLTPHOSRawAnalyzerCrudeComponent: public AliHLTPHOSRawAnalyzerComponentv3
{
 public:
  AliHLTPHOSRawAnalyzerCrudeComponent();
  virtual ~AliHLTPHOSRawAnalyzerCrudeComponent();
  AliHLTPHOSRawAnalyzerCrudeComponent(const AliHLTPHOSRawAnalyzerCrudeComponent & );
  AliHLTPHOSRawAnalyzerCrudeComponent & operator = (const AliHLTPHOSRawAnalyzerCrudeComponent)
  {
    return *this;
  };

  virtual int Deinit();
  virtual const char* GetComponentID();
  virtual AliHLTComponent* Spawn();
};

#endif
