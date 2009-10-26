//-*- Mode: C++ -*-
// $Id: AliHLTCALORawAnalyzerPeakFinder.h 29824 2008-11-10 13:43:55Z richterm $

#ifndef ALIHLTCALORAWANALYZERPEAKFINDER_H
#define ALIHLTCALORAWANALYZERPEAKFINDER_H
/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

//#include <Rtypes.h>
//#include "TObject.h"
#include "AliHLTCaloRawAnalyzer.h"
//#include "AliHLTCaloBase.h"
//class AliHLTCaloUtilities;

class AliHLTCaloRawAnalyzerPeakFinder : public AliHLTCaloRawAnalyzer
{
 public:
  AliHLTCaloRawAnalyzerPeakFinder();
  //  AliHLTCaloRawAnalyzerPeakFinder(const AliHLTCaloRawAnalyzerPeakFinder & );
  //  AliHLTCaloRawAnalyzerPeakFinder & operator = (const AliHLTCaloRawAnalyzerPeakFinder &)
  //    {
  //     return *this; 
  //    }

  virtual ~AliHLTCaloRawAnalyzerPeakFinder();


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
//  virtual void SetTVector(Double_t *tVectPtr =0, Int_t size = 0);
  void SetTVector(Double_t *tVectPtr =0, Int_t size = 0);


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
//  virtual void SetAVector(Double_t *aVectPtr =0, Int_t size =0);
 void SetAVector(Double_t *aVectPtr =0, Int_t size =0);

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
 AliHLTCaloRawAnalyzerPeakFinder(const AliHLTCaloRawAnalyzerPeakFinder & );
 AliHLTCaloRawAnalyzerPeakFinder & operator = (const AliHLTCaloRawAnalyzerPeakFinder &);

 Double_t   *fTVectorPtr;  //[1008]        /**<Peakfinder vector for TOF reconstruction*/
 Double_t   *fAVectorPtr;  //[1008]        /**<Peakfinder vector for Energy reconstruction*/  
 Int_t       fTVectorSize;
 Int_t       fAVectorSize;
 
 // AliHLTCaloUtilities *fUtilitiesPtr;


  //  ClassDef(AliHLTCaloRawAnalyzerPeakFinder, 2) 
  
    };

#endif
