//-*- Mode: C++ -*-
// $Id$


#ifndef ALIHLTPHOSRAWANALYZERLMS_H
#define ALIHLTPHOSRAWANALYZERLMS_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSRawAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

class AliHLTPHOSRawAnalyzerLMS : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerLMS();

  /**
   * Main constructor
   * @param dataPtr Data array for wich a subarray will be taken to perform the fit
   * @param fs the sampling frequency in entities of MHz. Needed in order to calculate physical time
   **/
  AliHLTPHOSRawAnalyzerLMS(double *dataPtr, double fs);
  AliHLTPHOSRawAnalyzerLMS(const AliHLTPHOSRawAnalyzerLMS & );
  AliHLTPHOSRawAnalyzerLMS & operator = (const AliHLTPHOSRawAnalyzerLMS)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSRawAnalyzerLMS();


  /**
   * Extraction of timing and energy using the LMS method.
   * The. The parameters "start" and "length" defines a sub array  of the data array
   * that will be used for the the fit. If start+length must not exeed the total length
   * of the Data array. "start" must be chosen as close as possible to t0.
   * The baseline must also be subtracted.
   * The length of "tVector" and "aVector" mus be equal to length.
   * "index + length" must not exeed the length of the data array set in the constructor.
   * @param start the start index of the subarray of the data array. 
   * @param length the number of samples to use starting from index 
   **/
 
  virtual void Evaluate(int start = 0, int lenght = 100);

 private:
  double   fMCovarPtrPtr[1008][1008]; /**<Covariance matrix of the measurements*/
  double   fPCovarPtrPtr[1008][1008];   /**<Covariance matrix of the estimated parameters*/
  ClassDef(AliHLTPHOSRawAnalyzerLMS, 2) 
  
    };

#endif
