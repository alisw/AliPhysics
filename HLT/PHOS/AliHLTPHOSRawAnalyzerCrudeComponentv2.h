#ifndef ALIHLTPHOSRAWANALYZERCRUDECOMPONENTV2_H
#define ALIHLTPHOSRAWANALYZERCRUDECOMPONENTV2_H

#include "AliHLTPHOSRawAnalyzerComponentv2.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 



class AliHLTPHOSRawAnalyzerCrudeComponentv2: public AliHLTPHOSRawAnalyzerComponentv2
{
 public:
  AliHLTPHOSRawAnalyzerCrudeComponentv2();
  virtual ~AliHLTPHOSRawAnalyzerCrudeComponentv2();
  AliHLTPHOSRawAnalyzerCrudeComponentv2(const AliHLTPHOSRawAnalyzerCrudeComponentv2 & );
  AliHLTPHOSRawAnalyzerCrudeComponentv2 & operator = (const AliHLTPHOSRawAnalyzerCrudeComponentv2)
  {
    return *this;
  };

  virtual int Deinit();
  virtual const char* GetComponentID();
  virtual AliHLTComponent* Spawn();
};

#endif
