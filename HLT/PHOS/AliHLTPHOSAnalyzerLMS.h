#ifndef ALIHLTPHOSANALYZERLMS_H
#define ALIHLTPHOSANALYZERLMS_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */


class AliHLTPHOSAnalyzerLMS : public AliHLTPHOSAnalyzer
{
 public:
  AliHLTPHOSAnalyzerLMS();
  AliHLTPHOSAnalyzerLMS(double *dataPtr, double fs);
  AliHLTPHOSAnalyzerLMS(const AliHLTPHOSAnalyzerLMS & );
  AliHLTPHOSAnalyzerLMS & operator = (const AliHLTPHOSAnalyzerLMS)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSAnalyzerLMS();
  virtual void Evaluate(int start = 0, int lenght = 100);
 private:
  double   **kfMCovarPtrPtr;   /**<Covariance matrix of the measurements*/
  double   **fPCovarPtrPtr;    /**<Covariance matrix of the estimated parameters*/

  ClassDef(AliHLTPHOSAnalyzerLMS, 2) 
  
    };

#endif
