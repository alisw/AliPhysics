#ifndef ALIPHOSPULSEGENERATOR_H
#define ALIPHOSPULSEGENERATOR_H
/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

class AliPHOSPulseGenerator
{
 public:
  /**
   * Default constructor
   **/
  AliPHOSPulseGenerator();
  
  /**
   * Destructor
   **/
  ~AliPHOSPulseGenerator();

  /**
   * Contruct a pulsegenrator object an initializes all necessary parameters
   * @param a Amplitude in ADC levels (0 -1023)
   * @param t0 Timedelay in nanoseconds of signal relative the first sample. This value should be between 0 and Ts
   * where Ts is the sample interval
   **/
  AliPHOSPulseGenerator(double a, double t0, const int N, const double t, const double f);

  /**
   * Adds a baseline offset to the signal
   * @param baselineLevel The basline level to add
   * @param *samples The sample array for which to add te basline offset
   **/
  void AddBaseline(double baselineLevel, double *samples);

  /**
   * Adds pretrigger samples to the sample array and returns 
   * a new array containing the pretrigger samples concatenatet
   * in front of the samples given by "samples"
   * @param The baseline value of the pretrigger samples
   * @param The sample array for which to add the pretrigger samples
   **/
  double *AddPretriggerSamples(double baslineLevel, double *samples);

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
   * @param -30DB cutoff frequency of the noise in entities of sampling frequency
   **/
  void AddNoise(double *dataPtr, double *sigma, double cutoff);

  /**
   * Returns the generated pulse with the parameters given in the constructor
   **/
  double *GetPulse();

  /**
   * Returns a Pulse with new amplidude and t0
   * @param a new amplidude, overriding the one given in the constructor
   **/
  double *GetPulse(double a, double t0);

  /**
   * Emulates the ADC. Rounds down to nearest Integerevalue all entries given by
   * dataPtr
   **/
  void Quantisize(double *dataPtr);

  /** sets the timedelay. Take as argument the timedelay in
   * nanoseconds and sets in i entities of sample intervals
   *   
   **/
 
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
  double *fDataPtr;
  double  fDT;
  double *fEvent;
};

#endif

