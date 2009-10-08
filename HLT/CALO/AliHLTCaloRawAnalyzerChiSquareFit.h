//-*- Mode: C++ -*-
// $Id: AliHLTCALORawAnalyzerChiSquareFit.h 30036 2008-11-24 16:43:44Z odjuvsla $

#ifndef ALIHLTCALORAWANALYZERCHISQUAREFIT_H
#define ALIHLTCALORAWANALYZERCHISQUAREFIT_H
//#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTCaloRawAnalyzer.h"

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/**
 * The AliHLTCaloPeakfinder class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the CALO raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/

class AliHLTCaloRawAnalyzerChiSquareFit : public AliHLTCaloRawAnalyzer
{
 public:

  AliHLTCaloRawAnalyzerChiSquareFit();
  AliHLTCaloRawAnalyzerChiSquareFit(const AliHLTCaloRawAnalyzerChiSquareFit & );

  AliHLTCaloRawAnalyzerChiSquareFit & operator = (const AliHLTCaloRawAnalyzerChiSquareFit)
    {
      return *this; 
    }
  
  virtual ~AliHLTCaloRawAnalyzerChiSquareFit();
  
  /**
   * Extraction of timing and energy using Chi Square Fit.
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

  //  ClassDef(AliHLTCaloRawAnalyzerChiSquareFit, 2) 
  
};

#endif
