#ifndef ALIHLTPHOSRAWANALYZERPEAKFINDER_H
#define ALIHLTPHOSRAWANALYZERPEAKFINDER_H
/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSRawAnalyzer.h"


class AliHLTPHOSRawAnalyzerPeakFinder : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerPeakFinder();
  AliHLTPHOSRawAnalyzerPeakFinder(const AliHLTPHOSRawAnalyzerPeakFinder & );
  AliHLTPHOSRawAnalyzerPeakFinder & operator = (const AliHLTPHOSRawAnalyzerPeakFinder)
    {
      return *this; 
    }

  
  virtual ~AliHLTPHOSRawAnalyzerPeakFinder();
  void SetTVector(double *tVector);
  void SetAVector(double *aVector);
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:
  double    *tVector;  //[1008]        /**<Peakfinder vector for TOF reconstruction*/
  double    *aVector;  //[1008]        /**<Peakfinder vector for Energy reconstruction*/

  ClassDef(AliHLTPHOSRawAnalyzerPeakFinder, 2) 
  
    };

#endif
