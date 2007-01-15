#ifndef ALIHLTPHOSPULSEGENERATOR_H
#define ALIHLTPHOSPULSEGENERATOR_H
/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

//#include <TObject.h>
#include <Rtypes.h>

class AliHLTPHOSPulseGenerator
{
 public:
  AliHLTPHOSPulseGenerator();
  virtual ~AliHLTPHOSPulseGenerator();
  AliHLTPHOSPulseGenerator(double a, double t0, const int N, const double t, const double f);
  AliHLTPHOSPulseGenerator(const AliHLTPHOSPulseGenerator & );
  AliHLTPHOSPulseGenerator & operator = (const AliHLTPHOSPulseGenerator)
    {
      return *this; 
    }
  void AddBaseline(double baselineLevel, double *samples);
  void AddNoise(double *dataPtr, double *sigma);
  void AddNoise(double *dataPtr, double *sigma, double cutoff);
  double *AddPretriggerSamples(double baslineLevel, double *samples);
  double *GetPulse();
  double *GetPulse(double a, double t0);
  void Quantisize(double *dataPtr);
  void SetAmplitude(double a);
  void SetDT(double fs);
  void SetTZero(double t0);

 private:
  void MakePulse(double *dtaPtr);
  void MakePulse(double *dtaPtr, double ampl);  
  double  fAmplitude;
  int     fNSamples;
  double  fTau;
  double  fSampleFreq;
  double  fTZero;
  double *fDataPtr; //[1000]
  double  fDT;
  double *fEvent; //[1000]
  
  ClassDef(AliHLTPHOSPulseGenerator,1)

};

#endif

