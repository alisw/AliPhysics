/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE Off-line Project.               *
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

#include "AliPHOSPulseGenerator.h"
//#include <stdio.h>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl; 

AliPHOSPulseGenerator::AliPHOSPulseGenerator()
{
  cout << "You cannot invoke the Pulsgenerator without parameters" << endl;
  fDataPtr=0;
}

AliPHOSPulseGenerator::~AliPHOSPulseGenerator()
{
  delete fDataPtr;
  fDataPtr=0;
}

AliPHOSPulseGenerator::AliPHOSPulseGenerator(double a, double t0, int N, double t, double fs)
{
  fDataPtr = new double[N];
  SetAmplitude(a);
  SetDT(fs);
  SetTZero(t0);
  fNSamples=N;
  fTau=t;
  fSampleFreq=fs;
  //  dT=tau/fs;   //Function values are calculated at intervals dT
  //  fDT=1/fs;   //Function values are calculated at intervals dT
  fDataPtr = new double[N];
  MakePulse(fDataPtr);
}

void 
AliPHOSPulseGenerator::AddBaseline(double baselineLevel, double *samples)
{
  cout << "AddBaseline not implemented yet" << endl;
}

void 
AliPHOSPulseGenerator::AddNoise(double *dataPtr, double *sigma)
{
  cout << "AddNoise is not implemented yet" << endl;
}

void 
AliPHOSPulseGenerator::AddNoise(double *dataPtr, double *sigma, double cutoff)
{
  cout << "AddNoise is not implemeted yet" << endl;
}


double *
AliPHOSPulseGenerator::AddPretriggerSamples(double baslineLevel, double *samples)
{
  cout << "AddPretriggerSamples not implemented yet" << endl;
}



double *
AliPHOSPulseGenerator::GetPulse()
{
  return fDataPtr;
}

void 
AliPHOSPulseGenerator::Quantisize(double *dataPtr)
{
  //  cout << "Quantisize is not implemented yet" << endl;
}

void
AliPHOSPulseGenerator::SetAmplitude(double a)
{
  fAmplitude=a;
}

void 
AliPHOSPulseGenerator::SetDT(double fs)
{
  fDT=1/fs;  
}

void
AliPHOSPulseGenerator::SetTZero(double t0)
{
  fTZero = -t0/1000; // Since time is in nanoseconds and the samplingfrequency is in MHz -> divide by 1000
}

//private
void
AliPHOSPulseGenerator::MakePulse(double *dtaPtr)
{
  
for(int i=0; i<fNSamples; i++)
  {
    dtaPtr[i]=fAmplitude*exp(2)*pow((i*fDT-fTZero)/fTau, 2)*exp(-2*(i*fDT-fTZero)/fTau);
  }  
}
