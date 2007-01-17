#ifndef ALIHLTPHOSANALYZERLMS_H
#define ALIHLTPHOSANALYZERLMS_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSRawAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


class AliHLTPHOSRawAnalyzerLMS : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerLMS();
  AliHLTPHOSRawAnalyzerLMS(double *dataPtr, double fs);
  AliHLTPHOSRawAnalyzerLMS(const AliHLTPHOSRawAnalyzerLMS & );
  AliHLTPHOSRawAnalyzerLMS & operator = (const AliHLTPHOSRawAnalyzerLMS)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSRawAnalyzerLMS();
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:
  double   **kfMCovarPtrPtr;   /**<Covariance matrix of the measurements*/
  double   **fPCovarPtrPtr;    /**<Covariance matrix of the estimated parameters*/

  ClassDef(AliHLTPHOSRawAnalyzerLMS, 2) 
  
    };

#endif
