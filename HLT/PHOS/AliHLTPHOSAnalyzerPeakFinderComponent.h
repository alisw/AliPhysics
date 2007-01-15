#ifndef ALIHLTPHOSANALYZERPEAKFINDERCOMPONENT_H
#define ALIHLTPHOSANALYZERPEAKFINDERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTPHOSAnalyzerComponent.h"

class AliHLTPHOSAnalyzerPeakFinderComponent: public AliHLTPHOSAnalyzerComponent
{
  AliHLTPHOSAnalyzerPeakFinderComponent();
  ~AliHLTPHOSAnalyzerPeakFinderComponent();
  AliHLTPHOSAnalyzerPeakFinderComponent(const AliHLTPHOSAnalyzerPeakFinderComponent & );
  AliHLTPHOSAnalyzerPeakFinderComponent & operator = (const AliHLTPHOSAnalyzerPeakFinderComponent)
  {
    return *this;
  };

ClassDef(AliHLTPHOSAnalyzerPeakFinderComponent, 2) 
  };



#endif
