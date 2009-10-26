//-*- Mode: C++ -*-
#ifndef ALIHLTCALORAWANALYZER_H
#define ALIHLTCALORAWANALYZER_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliHLTCALORawAnalyzer.h 34264 2009-08-14 18:29:23Z odjuvsla $ */

#include "Rtypes.h"
class AliHLTCaloUtilities;

class AliHLTCaloRawAnalyzer
{
 public:
  AliHLTCaloRawAnalyzer();
  virtual ~AliHLTCaloRawAnalyzer();
  void SetCorrectBaselineUsingFirstFiveSamples();
  void CorrectBaselineUsingFirstFiveSamples(UShort_t *data, const int length);
  inline float GetTiming() const {  return fDTof;};
  inline float GetEnergy() const { return fDAmpl;};
  void SetData(const UShort_t *data, const int length);
  virtual void Evaluate(Int_t start = 0, Int_t lenght = 100) = 0;

 protected:
  bool fDoCorrectBaselineUsingFirstFiveSamples;
  UShort_t   *fShortDataPtr;     /**<data that should be fitted */
  double     fSampleFrequency; /**<The ADC sample frequency in MHz used under data taking */
  double     fTau;	       /**<The risetime in micro seconds*/		 
  double     fDTof;            /**<Time of flight in entities of sample intervals */
  double     fDAmpl;           /**<Amplitude in entities of ADC levels*/

 protected:
  AliHLTCaloUtilities *fUtilitiesPtr;

 private:
  AliHLTCaloRawAnalyzer(const AliHLTCaloRawAnalyzer & );
  AliHLTCaloRawAnalyzer & operator = (const AliHLTCaloRawAnalyzer &);
};


#endif
