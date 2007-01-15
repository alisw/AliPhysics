#ifndef ALIHLTPHOSANALYZERCRUDECOMPONENT_H
#define ALIHLTPHOSANALYZERCRUDECOMPONENT_H

#include "AliHLTPHOSAnalyzerComponent.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 


class AliHLTPHOSAnalyzerCrudeComponent: public AliHLTPHOSAnalyzerComponent
{
  AliHLTPHOSAnalyzerCrudeComponent();
  ~AliHLTPHOSAnalyzerCrudeComponent();
  AliHLTPHOSAnalyzerCrudeComponent(const AliHLTPHOSAnalyzerCrudeComponent & );
  AliHLTPHOSAnalyzerCrudeComponent & operator = (const AliHLTPHOSAnalyzerCrudeComponent)
  {
    return *this;
  };

  //ClassDef(AliHLTPHOSAnalyzerCrudeComponent, 2) 
  };

#endif
