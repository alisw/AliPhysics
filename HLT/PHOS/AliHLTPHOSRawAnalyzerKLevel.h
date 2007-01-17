#ifndef ALIHLTPHOSRAWANALYZERKLEVEL_H
#define ALIHLTPHOSRAWANALYZERKLEVEL_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSRawAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

class AliHLTPHOSRawAnalyzerKLevel : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerKLevel();
  AliHLTPHOSRawAnalyzerKLevel(const AliHLTPHOSRawAnalyzerKLevel & );
  AliHLTPHOSRawAnalyzerKLevel & operator = (const AliHLTPHOSRawAnalyzerKLevel)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSRawAnalyzerKLevel();
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:
  double tKLevel;
  ClassDef(AliHLTPHOSRawAnalyzerKLevel, 2) 
  
    };

#endif
