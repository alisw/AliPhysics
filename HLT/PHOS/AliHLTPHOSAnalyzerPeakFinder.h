#ifndef ALIHLTPHOSANALYZERPEAKFINDER_H
#define ALIHLTPHOSANALYZERPEAKFINDER_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

class AliHLTPHOSAnalyzerPeakFinder : public AliHLTPHOSAnalyzer
{
 public:
  AliHLTPHOSAnalyzerPeakFinder();
  AliHLTPHOSAnalyzerPeakFinder(const AliHLTPHOSAnalyzerPeakFinder & );
  AliHLTPHOSAnalyzerPeakFinder & operator = (const AliHLTPHOSAnalyzerPeakFinder)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSAnalyzerPeakFinder();
  void SetTVector(double *tVector);
  void SetAVector(double *aVector);
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:
  double    *tVector;          /**<Peakfinder vector for TOF reconstruction*/
  double    *aVector;          /**<Peakfinder vector for Energy reconstruction*/

  ClassDef(AliHLTPHOSAnalyzerPeakFinder, 2) 
  
    };

#endif
