//-*- Mode: C++ -*-
#ifndef ALIHLTCALORAWANALYZER_H
#define ALIHLTCALORAWANALYZER_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliHLTCALORawAnalyzer.h 34264 2009-08-14 18:29:23Z odjuvsla $ */

#include "Rtypes.h"


//#include "AliHLTCaloBase.h"

class AliHLTCaloUtilities;

class AliHLTCaloRawAnalyzer
//: public AliHLTCaloBase
//class AliHLTCaloRawAnalyzer
{
 public:
  AliHLTCaloRawAnalyzer();
  virtual ~AliHLTCaloRawAnalyzer();
  AliHLTCaloRawAnalyzer(double *dataPtr, double fs);

  void SetCorrectBaselineUsingFirstFiveSamples();
  void CorrectBaselineUsingFirstFiveSamples(UInt_t *data, const int length);
  // void CorrectBaselineUsingFirstFiveSamples(int *data, int length);
  //void CorrectBaselineUsingFirstFiveSamples(double *data, int length);


  //  void BaselineCorrection(double *dataPtr, int N);
  void BaselineCorrection(double *dataPtr, double baselineValue);  

  int FindStartIndex(double treshold);
  float GetTiming() const;
  float GetEnergy() const;

  //  void SetData(const double *data);
  void SetData(const UInt_t *data, const int length);
  void SetData(const UShort_t *data, const int length);
  // void SetData(UInt_t *data, const int length);

  // void SetData(double *data, const int length);
  //  void SetData(int *data, const int length);

 
  void SetSampleFreq(double freq);
  void SetStartIndex(int startIndex);
  void MakeInitialGuess();
  void MakeInitialGuess(int treshold);
  
  virtual void SetTVector(Double_t *tVector, Int_t size);
  virtual void SetAVector(Double_t *aVector, Int_t size);

  virtual void Evaluate(Int_t start = 0, Int_t lenght = 100) = 0;


 protected:
  bool fDoCorrectBaselineUsingFirstFiveSamples;

  //  mutable  Double_t   *fFloatDataPtr;   /**<Float representation of data that should be fitted */
  //  mutable  UInt_t     *fIntDataPtr;     /**<data that should be fitted */

 /*  mutable const Double_t   *fFloatDataPtr;   /\**<Float representation of data that should be fitted *\/ */
/*   mutable const unsigned int     *fIntDataPtr;     /\**<data that should be fitted *\/ */


  //double   *fFloatDataPtr;   /**<Float representation of data that should be fitted */

  double   *fDoubleDataPtr;   /**<Float representation of data that should be fitted */
  UInt_t   *fIntDataPtr;     /**<data that should be fitted */
  UShort_t   *fShortDataPtr;     /**<data that should be fitted */
 //  int      *fIntDataPtr;     /**<data that should be fitted */


  //  mutable const UInt_t     *fIntDataPtr;     /**<data that should be fitted */
  
  double     fSampleFrequency; /**<The ADC sample frequency in MHz used under data taking */
  double     fDTofGuess;       /**<Initial guess for t0*/
  double     fDAmplGuess;      /**<Initial guess for amplitude*/
  double     fTau;	       /**<The risetime in micro seconds*/		 
  double     fDTof;            /**<Time of flight in entities of sample intervals */
  double     fDAmpl;           /**<Amplitude in entities of ADC levels*/
  int        fStartIndex;      /**<Starindex of the time dependent altro signal*/

  //  double unsigned fTest;

  bool fUseShortValues;

 protected:
  AliHLTCaloUtilities *fUtilitiesPtr;

 private:
  AliHLTCaloRawAnalyzer(const AliHLTCaloRawAnalyzer & );
  AliHLTCaloRawAnalyzer & operator = (const AliHLTCaloRawAnalyzer &);
};


#endif
