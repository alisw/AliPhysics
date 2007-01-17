#ifndef ALIHLTPHOSANALYZERCRUDECOMPONENT_H
#define ALIHLTPHOSANALYZERCRUDECOMPONENT_H

#include "AliHLTPHOSRawAnalyzerComponent.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 


class AliHLTPHOSRawAnalyzerCrudeComponent: public AliHLTPHOSRawAnalyzerComponent
{
  AliHLTPHOSRawAnalyzerCrudeComponent();
  ~AliHLTPHOSRawAnalyzerCrudeComponent();
  AliHLTPHOSRawAnalyzerCrudeComponent(const AliHLTPHOSRawAnalyzerCrudeComponent & );
  AliHLTPHOSRawAnalyzerCrudeComponent & operator = (const AliHLTPHOSRawAnalyzerCrudeComponent)
  {
    return *this;
  };

  //ClassDef(AliHLTPHOSRawAnalyzerCrudeComponent, 2) 
  };

#endif
