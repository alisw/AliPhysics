#ifndef ALIHLTPHOSANALYZERCRUDE_H
#define ALIHLTPHOSANALYZERCRUDE_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSAnalyzer.h"


/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


class AliHLTPHOSAnalyzerCrude : public AliHLTPHOSAnalyzer
{
 public:
  AliHLTPHOSAnalyzerCrude();
  AliHLTPHOSAnalyzerCrude(const AliHLTPHOSAnalyzerCrude & );
  AliHLTPHOSAnalyzerCrude & operator = (const AliHLTPHOSAnalyzerCrude)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSAnalyzerCrude();
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:
  ClassDef(AliHLTPHOSAnalyzerCrude, 2) 
  
    };

#endif
