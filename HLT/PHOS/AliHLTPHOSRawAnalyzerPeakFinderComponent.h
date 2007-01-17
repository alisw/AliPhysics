#ifndef ALIHLTPHOSANALYZERPEAKFINDERCOMPONENT_H
#define ALIHLTPHOSANALYZERPEAKFINDERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTPHOSRawAnalyzerComponent.h"

class AliHLTPHOSRawAnalyzerPeakFinderComponent: public AliHLTPHOSRawAnalyzerComponent
{
  AliHLTPHOSRawAnalyzerPeakFinderComponent();
  ~AliHLTPHOSRawAnalyzerPeakFinderComponent();
  AliHLTPHOSRawAnalyzerPeakFinderComponent(const AliHLTPHOSRawAnalyzerPeakFinderComponent & );
  AliHLTPHOSRawAnalyzerPeakFinderComponent & operator = (const AliHLTPHOSRawAnalyzerPeakFinderComponent)
  {
    return *this;
  };

  //ClassDef(AliHLTPHOSRawAnalyzerPeakFinderComponent, 2) 
  };



#endif
