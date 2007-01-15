#ifndef ALIHLTPHOSANALYZERKLEVEL_H
#define ALIHLTPHOSANALYZERKLEVEL_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

class AliHLTPHOSAnalyzerKLevel : public AliHLTPHOSAnalyzer
{
 public:
  AliHLTPHOSAnalyzerKLevel();
  AliHLTPHOSAnalyzerKLevel(const AliHLTPHOSAnalyzerKLevel & );
  AliHLTPHOSAnalyzerKLevel & operator = (const AliHLTPHOSAnalyzerKLevel)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSAnalyzerKLevel();
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:
  double tKLevel;
  ClassDef(AliHLTPHOSAnalyzerKLevel, 2) 
  
    };

#endif
