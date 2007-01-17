#ifndef ALIHLTPHOSRAWANALYZERCHISQUAREFIT_H
#define ALIHLTPHOSRAWANALYZERCHISQUAREFIT_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSRawAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


class AliHLTPHOSRawAnalyzerChiSquareFit : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerChiSquareFit();
  AliHLTPHOSRawAnalyzerChiSquareFit(const AliHLTPHOSRawAnalyzerChiSquareFit & );

  AliHLTPHOSRawAnalyzerChiSquareFit & operator = (const AliHLTPHOSRawAnalyzerChiSquareFit)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSRawAnalyzerChiSquareFit();
  virtual void Evaluate(int start = 0, int lenght = 100);
   

 private:

  ClassDef(AliHLTPHOSRawAnalyzerChiSquareFit, 2) 
  
    };

#endif
