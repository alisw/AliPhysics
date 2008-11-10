//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSPULSEGENERATOR_H
#define ALIHLTPHOSPULSEGENERATOR_H
/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include <Rtypes.h>


class AliHLTPHOSPulseGenerator
{
 public:
  AliHLTPHOSPulseGenerator();
  virtual ~AliHLTPHOSPulseGenerator();

  /**
   * Contruct a pulsegenrator object an initializes all necessary parameters
   * @param a Amplitude in ADC levels (0 -1023)
   * @param t0 Timedelay in nanoseconds of signal relative the first sample. This value should be between 0 and Ts
   * @param N the number of samples
   * @param tau Rise time of the semi Gaussian signal
   * @param fs samling rate
   **/
  AliHLTPHOSPulseGenerator(double a, double t0, const int N , const double tau, const double fs);

  AliHLTPHOSPulseGenerator(const AliHLTPHOSPulseGenerator & );


  AliHLTPHOSPulseGenerator & operator = (const AliHLTPHOSPulseGenerator &)
    {
      return *this; 
    }

  /**
   * Adds a baseline offset to the signal
   * @param baselineLevel The basline level to add
   * @param *samples The sample array for which to add te basline offset
   **/
  void AddBaseline(double baselineLevel = 0, double *samples = 0);

  /**
   * Adds Gaussian white noise to the sample array given by *dataPtr.
   * @param dataPtr array of samples
   * @param sigma the noise amplitude in entities of ADC levels  
   **/
  void AddNoise(double *dataPtr, double *sigma);

  /**
   * Adds correlated Gaussian noise with cutof frequency "cutoff"
   * @param dataPtr array of values
   * @param sigma noise amplitude in entities of ADC levels
   * @param cutoff -30DB cutoff frequency of the noise in entities of sampling frequency
   **/
  void AddNoise(double *dataPtr, double *sigma, double cutoff); 


  /**
   * Adds pretrigger samples to the sample array and returns 
   * a new array containing the pretrigger samples concatenatet
   * in front of the samples given by "samples"
   * @param baselineLevel The baseline value of the pretrigger samples
   * @param samples The sample array for which to add the pretrigger samples
   **/
  double *AddPretriggerSamples(double baselineLevel = 0, double *samples = 0);


  /**
   * Returns a Pulse with new amplidude and t0
   * @param a new amplidude, overriding the one given in the constructor
   * @param t0 start time of the pulse relative to the sampling clock.
   **/
  double *GetPulse(double a = 1, double t0 = 0);


  /**
   * Emulates the ADC. Rounds down to nearest Integerevalue all entries given by
   * dataPtr
   **/
  void Quantisize(double *dataPtr) const;

  void SetAmplitude(double a = 1);
  void SetDT(double fs = 10);
  void SetTZero(double t0 = 0);
  void SetSampleFreq(int fs);
  void MakePulse(double *dtaPtr, int N);
  void MakePulse(double *dtaPtr);

 private:
  double  fAmplitude;         /**<The amplitude in entities of ADC counts of the genrated pulse*/
  int     fNSamples;          /**<The number of samples of the genrated pulse*/
  double  fTau;               /**<The risetime in entities of us of the generated pulse*/
  double  fSampleFreq;        /**<The sampling frequency in MHz*/ 
  double  fTZero;             /**<t0 of the genrated pulse in entities of nanoseconds*/   
  double *fDataPtr; //[1000]  /**<pointer to array holding the genrated pulse*/ 
  double  fDT;                /**<1/fSampleFreq*/
  //  double *fEvent; //[1000]
  ClassDef(AliHLTPHOSPulseGenerator,1)
};

#endif

