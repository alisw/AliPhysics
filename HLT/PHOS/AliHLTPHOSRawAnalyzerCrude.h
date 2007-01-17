#ifndef ALIHLTPHOSRAWANALYZERCRUDE_H
#define ALIHLTPHOSRAWANALYZERCRUDE_H
#include <Rtypes.h>
#include "TObject.h"
//         "AliHLTPHOSRawAnalyzer"
#include "AliHLTPHOSRawAnalyzer.h"


/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


class AliHLTPHOSRawAnalyzerCrude : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerCrude();
  AliHLTPHOSRawAnalyzerCrude(const AliHLTPHOSRawAnalyzerCrude & );
  AliHLTPHOSRawAnalyzerCrude & operator = (const AliHLTPHOSRawAnalyzerCrude)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSRawAnalyzerCrude();
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:
  ClassDef(AliHLTPHOSRawAnalyzerCrude, 2) 
  
    };

#endif
