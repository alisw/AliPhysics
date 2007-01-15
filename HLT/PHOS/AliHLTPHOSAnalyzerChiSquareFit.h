#ifndef ALIHLTPHOSANALYZERCHISQUAREFIT_H
#define ALIHLTPHOSANALYZERCHISQUAREFIT_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


class AliHLTPHOSAnalyzerChiSquareFit : public AliHLTPHOSAnalyzer
{
 public:
  AliHLTPHOSAnalyzerChiSquareFit();
  AliHLTPHOSAnalyzerChiSquareFit(const AliHLTPHOSAnalyzerChiSquareFit & );

  AliHLTPHOSAnalyzerChiSquareFit & operator = (const AliHLTPHOSAnalyzerChiSquareFit)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSAnalyzerChiSquareFit();
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:

  ClassDef(AliHLTPHOSAnalyzerChiSquareFit, 2) 
  
    };

#endif
