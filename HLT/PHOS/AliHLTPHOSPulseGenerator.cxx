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



AliHLTPHOSPulseGenerator::AliHLTPHOSPulseGenerator(): fAmplitude(0), fNSamples(0),fTau(0), fSampleFreq(0), fTZero(0), fDataPtr(0), fDT(0)
{
  //Never to be called   
  cout << "You cannot invoke the Pulsgenerator without parameters" << endl;
}


AliHLTPHOSPulseGenerator::AliHLTPHOSPulseGenerator(const AliHLTPHOSPulseGenerator &): fAmplitude(0), fNSamples(0),fTau(0), fSampleFreq(0), fTZero(0), fDataPtr(0), fDT(0)
{
  
}


AliHLTPHOSPulseGenerator::~AliHLTPHOSPulseGenerator()
{
  //Destructor  
  delete fDataPtr;
  fDataPtr=0;
}


AliHLTPHOSPulseGenerator::AliHLTPHOSPulseGenerator(double a, double t0, int N, double tau, double fs): fAmplitude(a), fNSamples(N),fTau(0), fSampleFreq(fs), fTZero(0), fDataPtr(0), fDT(0)
{
  //See header file for documentation  
  fDataPtr = new double[100];
  SetAmplitude(a);
  SetDT(fs);
  SetTZero(t0);
  fNSamples=N;
  fTau=tau;
  fSampleFreq=fs;
}


void 
AliHLTPHOSPulseGenerator::AddBaseline(double baselineLevel, double *samples)
{
  //See header file for documentation
  double *tmpSamples;
  tmpSamples = samples;
  printf("\nbaselineLevel = %f\n", baselineLevel);
  cout << "AddBaseline not implemented yet" << endl;
}


void 
AliHLTPHOSPulseGenerator::AddNoise(double *dataPtr, double *sigma)
{
  //See header file for documentation 
  printf("\ndataPtr = %f, sigma = %f\n", *dataPtr, *sigma);
  cout << "AddNoise is not implemented yet" << endl;
}


void 
AliHLTPHOSPulseGenerator::AddNoise(double *dataPtr, double *sigma, double cutoff)
{
  //See header file for documentation 
  printf("\ndataPtr = %f, sigma = %f, cutoff = %f\n", *dataPtr, *sigma, cutoff);
  cout << "AddNoise is not implemeted yet" << endl;
}


double *
AliHLTPHOSPulseGenerator::AddPretriggerSamples(double baselineLevel, double *samples)
{
    //See header file for documentation
  printf("\nbaslinelevel = %f, samples = %f\n", baselineLevel, *samples);
  cout << "AddPretriggerSamples not implemented yet" << endl;
  return 0;
}


double *
AliHLTPHOSPulseGenerator::GetPulse(double /*a*/, double /*t0*/)
{
  //See header file for documentation
  return fDataPtr;
}


void 
AliHLTPHOSPulseGenerator::Quantisize(double *dataPtr) const
{
  //See header file for documentation
  double *dtaPtr;
  dtaPtr = new double[100];
  dtaPtr = dataPtr;
}


void
AliHLTPHOSPulseGenerator::SetAmplitude(double a)
{
  //See header file for documentation
  fAmplitude=a;
}


void 
AliHLTPHOSPulseGenerator::SetDT(double fs)
{
  //See header file for documentation
  fDT=1/fs;  
}


void
AliHLTPHOSPulseGenerator::SetTZero(double t0)
{
  //See header file for documentation
  fTZero = -t0/1000; // Since time is in nanoseconds and the samplingfrequency is in MHz -> divide by 1000
}


void 
AliHLTPHOSPulseGenerator::SetSampleFreq(int fs)
{
  //See header file for documentation
  fSampleFreq = fs;
  SetDT(fs);
}


void
AliHLTPHOSPulseGenerator::MakePulse(double *dtaPtr)
{
  //See header file for documentation
  for(int i=0; i<fNSamples; i++)
    {
      dtaPtr[i]=fAmplitude*exp((Double_t)2)*pow((i*fDT-fTZero)/fTau, 2)*exp(-2*(i*fDT-fTZero)/fTau);
    }  
}


void
AliHLTPHOSPulseGenerator::MakePulse(double *dtaPtr, int N)
{
  //See header file for documentation
  for(int i=0; i<N; i++)
    {
      dtaPtr[i]=fAmplitude*exp((Double_t)2)*pow((i*fDT-fTZero)/fTau, 2)*exp(-2*(i*fDT-fTZero)/fTau);
    }  
}

