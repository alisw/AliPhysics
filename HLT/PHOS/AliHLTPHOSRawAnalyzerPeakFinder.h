//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRAWANALYZERPEAKFINDER_H
#define ALIHLTPHOSRAWANALYZERPEAKFINDER_H
/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliHLTPHOSRawAnalyzer.h"


class AliHLTPHOSRawAnalyzerPeakFinder : public AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzerPeakFinder();
  virtual ~AliHLTPHOSRawAnalyzerPeakFinder();


/**
* Extraction of timing and energy using the Peakfinde Algorithm.
* The. The parameters "start" and "length" defines a sub array  of the data array
* that will be used for the the fit. If start+length must not exeed the total length
* of the Data array. "start" must be chosen as close as possible to t0.
* The baseline must also be subtracted.
* The length of "tVector" and "aVector" mus be equal to length.
* "index + length" must not exeed the length of the data array set in the constructor.
* @param tVectPtr the peakfinder vector for timing
* @param size size in number of values of the time vector
*/
  virtual void SetTVector(Double_t *tVectPtr =0, Int_t size = 0);



/**
* Extraction of timing and energy using the Peakfinde Algorithm.
* The. The parameters "start" and "length" defines a sub array  of the data array
* that will be used for the the fit. If start+length must not exeed the total length
* of the Data array. "start" must be chosen as close as possible to t0.
* The baseline must also be subtracted.
* The length of "tVector" and "aVector" mus be equal to length.
* "index + length" must not exeed the length of the data array set in the constructor.
* @param aVectPtr the peakfinder vector for timing
* @param size size in number of values of the time vector
*/
  virtual void SetAVector(Double_t *aVectPtr =0, Int_t size =0);


/**
* Extraction of timing and energy using the Peakfinde Algorithm.
* The. The parameters "start" and "length" defines a sub array  of the data array
* that will be used for the the fit. If start+length must not exeed the total length
* of the Data array. "start" must be chosen as close as possible to t0.
* The baseline must also be subtracted.
* The length of "tVector" and "aVector" mus be equal to length.
* "index + length" must not exeed the length of the data array set in the constructor.
* @param start the start index of the subarray of the data array. 
* @param length the number of samples to use starting from index 
**/
  virtual void Evaluate(Int_t start = 0, Int_t length = 100);
 
 private:
 AliHLTPHOSRawAnalyzerPeakFinder(const AliHLTPHOSRawAnalyzerPeakFinder & );
 AliHLTPHOSRawAnalyzerPeakFinder & operator = (const AliHLTPHOSRawAnalyzerPeakFinder &);

 Double_t   *fTVectorPtr;  //[1008]        /**<Peakfinder vector for TOF reconstruction*/
 Double_t   *fAVectorPtr;  //[1008]        /**<Peakfinder vector for Energy reconstruction*/  
 Int_t       fTVectorSize;
 Int_t       fAVectorSize;
 
 // AliHLTPHOSUtilities *fUtilitiesPtr;


  ClassDef(AliHLTPHOSRawAnalyzerPeakFinder, 2) 
  
    };

#endif
