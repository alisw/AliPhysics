//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRAWANALYZERKLEVEL_H
#define ALIHLTPHOSRAWANALYZERKLEVEL_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSRawAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/**
 * The AliHLTPHOSPeakfinder class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the PHOS raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/

class AliHLTPHOSRawAnalyzerKLevel : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerKLevel();
  AliHLTPHOSRawAnalyzerKLevel(const AliHLTPHOSRawAnalyzerKLevel & );
  AliHLTPHOSRawAnalyzerKLevel & operator = (const AliHLTPHOSRawAnalyzerKLevel)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSRawAnalyzerKLevel();
 
  /**
   * Extraction of timing and energy using the K-Level method.
   * The. The parameters "start" and "length" defines a sub array  of the data array
   * that will be used for the the fit. If start+length must not exeed the total length
   * of the Data array. "start" must be chosen as close as possible to t0.
   * The baseline must also be subtracted.
   * The length of "tVector" and "aVector" mus be equal to length.
   * "index + length" must not exeed the length of the data array set in the constructor.
   * @param start the start index of the subarray of the data array. 
   * @param length the number of samples to use starting from index 
   **/
  virtual void Evaluate(int start = 0, int length = 100);
 private:
  double fTKLevel; /**<K-Level*/
  ClassDef(AliHLTPHOSRawAnalyzerKLevel, 2) 
  
    };

#endif
