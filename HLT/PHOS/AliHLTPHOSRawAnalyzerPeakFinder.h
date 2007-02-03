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
  virtual void SetTVector(Double_t *tVect, Int_t size);
  virtual void SetAVector(Double_t *aVect, Int_t size);
  virtual void Evaluate(Int_t start = 0, Int_t lenght = 100);
 private:
  Double_t   *fTVectorPtr;  //[1008]        /**<Peakfinder vector for TOF reconstruction*/
  Double_t   *fAVectorPtr;  //[1008]        /**<Peakfinder vector for Energy reconstruction*/  
  Int_t       fTVectorSize;
  Int_t       fAVectorSize;

  ClassDef(AliHLTPHOSRawAnalyzerPeakFinder, 2) 
  
    };

#endif
