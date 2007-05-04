/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSPulseGenerator.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::endl; 

ClassImp(AliHLTPHOSPulseGenerator) 

//______________________________________________________________________________________________________
AliHLTPHOSPulseGenerator::AliHLTPHOSPulseGenerator(): fAmplitude(0), fNSamples(0),fTau(0), fSampleFreq(0), fTZero(0), fDataPtr(0), fDT(0)
{
  cout << "You cannot invoke the Pulsgenerator without parameters" << endl;
}


//______________________________________________________________________________________________________
AliHLTPHOSPulseGenerator::AliHLTPHOSPulseGenerator(const AliHLTPHOSPulseGenerator &): fAmplitude(0), fNSamples(0),fTau(0), fSampleFreq(0), fTZero(0), fDataPtr(0), fDT(0)
{
  
}


//______________________________________________________________________________________________________
AliHLTPHOSPulseGenerator::~AliHLTPHOSPulseGenerator()
{
  delete fDataPtr;
  fDataPtr=0;
}

/**
 * Contruct a pulsegenrator object an initializes all necessary parameters
 * @param a Amplitude in ADC levels (0 -1023)
 * @param t0 Timedelay in nanoseconds of signal relative the first sample. This value should be between 0 and Ts
 * @param N the number of samples
 * @param tau Rise time of the semi Gaussian signal
 * @param fs samling rate
 **/
AliHLTPHOSPulseGenerator::AliHLTPHOSPulseGenerator(double a, double t0, int N, double tau, double fs): fAmplitude(a), fNSamples(N),fTau(0), fSampleFreq(fs), fTZero(0), fDataPtr(0), fDT(0)
{
  fDataPtr = new double[100];
  SetAmplitude(a);
  SetDT(fs);
  SetTZero(t0);
  fNSamples=N;
  fTau=tau;
  fSampleFreq=fs;
  //  MakePulse(fDataPtr,a);
  MakePulse(fDataPtr);
}


/**
 * Adds a baseline offset to the signal
 * @param baselineLevel The basline level to add
 * @param *samples The sample array for which to add te basline offset
 **/
void 
AliHLTPHOSPulseGenerator::AddBaseline(double baselineLevel, double *samples)
{
  double *tmpSamples;
  tmpSamples = samples;
  printf("\nbaselineLevel = %f\n", baselineLevel);
  cout << "AddBaseline not implemented yet" << endl;
}

/**
 * Adds Gaussian white noise to the sample array given by *dataPtr.
 * @param dataPtr array of samples
 * @param sigma the noise amplitude in entities of ADC levels  
 **/
void 
AliHLTPHOSPulseGenerator::AddNoise(double *dataPtr, double *sigma)
{
  printf("\ndataPtr = %f, sigma = %f\n", *dataPtr, *sigma);
  cout << "AddNoise is not implemented yet" << endl;
}


/**
 * Adds correlated Gaussian noise with cutof frequency "cutoff"
 * @param dataPtr array of values
 * @param sigma noise amplitude in entities of ADC levels
 * @param cutoff -30DB cutoff frequency of the noise in entities of sampling frequency
 **/
void 
AliHLTPHOSPulseGenerator::AddNoise(double *dataPtr, double *sigma, double cutoff)
{
  printf("\ndataPtr = %f, sigma = %f, cutoff = %f\n", *dataPtr, *sigma, cutoff);
  cout << "AddNoise is not implemeted yet" << endl;
}



/**
 * Adds pretrigger samples to the sample array and returns 
 * a new array containing the pretrigger samples concatenatet
 * in front of the samples given by "samples"
 * @param baselineLevel The baseline value of the pretrigger samples
 * @param samples The sample array for which to add the pretrigger samples
 **/
double *
AliHLTPHOSPulseGenerator::AddPretriggerSamples(double baselineLevel, double *samples)
{
  printf("\nbaslinelevel = %f, samples = %f\n", baselineLevel, *samples);
  cout << "AddPretriggerSamples not implemented yet" << endl;
  return 0;
}


/**
 * Returns a Pulse with new amplidude and t0
 * @param a new amplidude, overriding the one given in the constructor
 * @param t0 start time of the pulse relative to the sampling clock.
 **/
double *
AliHLTPHOSPulseGenerator::GetPulse(double a, double t0)
{
  return fDataPtr;
}

/**
 * Emulates the ADC. Rounds down to nearest Integerevalue all entries given by
 * dataPtr
 **/
void 
AliHLTPHOSPulseGenerator::Quantisize(double *dataPtr)
{
  double *dtaPtr;
  dtaPtr = new double[100];
  dtaPtr = dataPtr;
  //  cout << "Quantisize is not implemented yet" << endl;
}


//______________________________________________________________________________________________________
void
AliHLTPHOSPulseGenerator::SetAmplitude(double a)
{
  fAmplitude=a;
}


//______________________________________________________________________________________________________
void 
AliHLTPHOSPulseGenerator::SetDT(double fs)
{
  fDT=1/fs;  
}

//______________________________________________________________________________________________________
void
AliHLTPHOSPulseGenerator::SetTZero(double t0)
{
  fTZero = -t0/1000; // Since time is in nanoseconds and the samplingfrequency is in MHz -> divide by 1000
}


//______________________________________________________________________________________________________
void
AliHLTPHOSPulseGenerator::MakePulse(double *dtaPtr)
{
for(int i=0; i<fNSamples; i++)
  {
    dtaPtr[i]=fAmplitude*exp((Double_t)2)*pow((i*fDT-fTZero)/fTau, 2)*exp(-2*(i*fDT-fTZero)/fTau);
  }  
}
