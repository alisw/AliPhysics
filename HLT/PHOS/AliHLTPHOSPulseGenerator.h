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
  AliHLTPHOSPulseGenerator(double a, double t0, const int N , const double tau, const double fs);
  AliHLTPHOSPulseGenerator(const AliHLTPHOSPulseGenerator & );
  AliHLTPHOSPulseGenerator & operator = (const AliHLTPHOSPulseGenerator)
    {
      return *this; 
    }
  void AddBaseline(double baselineLevel = 0, double *samples = 0);
  void AddNoise(double *dataPtr, double *sigma);
  void AddNoise(double *dataPtr, double *sigma, double cutoff); 
  double *AddPretriggerSamples(double baslineLevel = 0, double *samples = 0);
  double *GetPulse(double a = 1, double t0 = 0);
  void Quantisize(double *dataPtr = 0);
  void SetAmplitude(double a = 1);
  void SetDT(double fs = 10);
  void SetTZero(double t0 = 0);

 private:
  void MakePulse(double *dtaPtr = 0);
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

