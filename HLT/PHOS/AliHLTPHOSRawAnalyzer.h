//-*- Mode: C++ -*-
#ifndef ALIHLTPHOSRAWANALYZER_H
#define ALIHLTPHOSRAWANALYZER_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

#include "Rtypes.h"

class AliHLTPHOSRawAnalyzer
{
 public:
  AliHLTPHOSRawAnalyzer();
  virtual ~AliHLTPHOSRawAnalyzer();
  int FindStartIndex(double treshold);
  float GetTiming() const { return fDTof;};  // peak position in entities of sample indexes
  float GetEnergy() const { return fDAmpl;};  // amplitude in entities of ADC channels   ; 
  void SetData(const UShort_t *data, const int /*length*/) {fShortDataPtr = const_cast<UShort_t *>(data);};
  void SetStartIndex(int startIndex) {fStartIndex = startIndex;};
  virtual void Evaluate(Int_t start = 0, Int_t lenght = 100) = 0;

protected:
  double   *fDoubleDataPtr;   /**<Float representation of data that should be fitted */
  UShort_t   *fShortDataPtr;     /**<data that should be fitted */
  double     fSampleFrequency; /**<The ADC sample frequency in MHz used under data taking */
  double     fDTofGuess;       /**<Initial guess for t0*/
  double     fDAmplGuess;      /**<Initial guess for amplitude*/
  double     fTau;	       /**<The risetime in micro seconds*/		 
  double     fDTof;            /**<Time of flight in entities of sample intervals */
  double     fDAmpl;           /**<Amplitude in entities of ADC levels*/
  int        fStartIndex;      /**<Starindex of the time dependent altro signal*/

private:
  AliHLTPHOSRawAnalyzer(const AliHLTPHOSRawAnalyzer & );
  AliHLTPHOSRawAnalyzer & operator = (const AliHLTPHOSRawAnalyzer &);
};


#endif
