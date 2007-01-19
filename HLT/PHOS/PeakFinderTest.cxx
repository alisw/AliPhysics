#include "AliHLTPHOSPulseGenerator.h"
#include "AliHLTPHOSRawAnalyzerPeakFinder.h"
#include <stdio.h>
#include <cmath>

void setFileName(char *fName, int start, int length, double tau, double fs);

/**
 * Testing of the Class AliPHOSFitter
 **/
int main()
{
  FILE *fp = 0;
  int start  = 0;               // Start index of subarray of sample array
  int N = 128;                  // Number of samples 
  char fileName[100];           // Name of file containing Peakfinder vectors
  double amplitude;             // Amplitude/energy in ADC levels
  double t0;                    // timedelay in nanoseconds  

  double tau = 2;               // risetime in microseconds
  double fs = 20;               // sample frequency in Megahertz
  double timeVector[N];         // Peakfinder vector for reconstruction of time
  double amplitudeVector[N];    // Peakfinder vector for reconstruction of energy
  double aSystError;
  double tSystError;

  int ts = (int)(1000/fs);  
  printf("\nts=%d\n", ts);

  printf("type amplitude in ADC levels (0-1023):");
  scanf("%lf", &amplitude);
  printf("type timedelay in nanoseconds (0-%d):", ts); 
  scanf("%lf", &t0);

  AliHLTPHOSPulseGenerator *pulseGenPtr = new AliHLTPHOSPulseGenerator(amplitude, t0, N, tau, fs);
  double *data = pulseGenPtr->GetPulse(); 

  setFileName(fileName, start, N, tau, fs);
  fp = fopen(fileName,"r");

  if(fp == 0)
    {
      printf("\nFile does not exist\n");
    }
  else
    {
      for(int i=0; i < N; i++)
	{
	  fscanf(fp, "%lf", &amplitudeVector[i]);
	}

      fscanf(fp, "\n");

     for(int i=0; i < N; i++)
	{
	  fscanf(fp, "%lf", &timeVector[i]);
	}

     fscanf(fp, "%lf", &aSystError);
     fscanf(fp, "%lf", &tSystError);
     printf("\nPeakfinder vectors loaded from %s\n", fileName);
    }


   tSystError = tSystError*pow(10, 9); //to give systematic error of timing in nanoseconds
   aSystError = aSystError*100;        //to give systematic error of amplitude in percent


   //  AliHLTPHOSAnalyzerPeakFinder *fitPtr= new AliHLTPHOSAnalyzerPeakFinder(data, fs); 
   AliHLTPHOSRawAnalyzerPeakFinder *fitPtr= new AliHLTPHOSRawAnalyzerPeakFinder(); 
  
   fitPtr->SetData(data);
   fitPtr->SetSampleFreq(fs);
   fitPtr->SetTVector(timeVector);
   fitPtr->SetAVector(amplitudeVector);
   //   fitPtr->Set
   fitPtr->Evaluate(start, N);
  
  double energy;
  double time;

  time    = fitPtr->GetTiming();
  energy  = fitPtr->GetEnergy();

  printf("\nReal amplitude \t\t= %lf ADC counts \nReconstructed amplitude\t= %lf ADC counts\n", amplitude, energy);
  printf("\nReal time \t\t= %lf nanoseconds \nReconstructed time\t= %lf nanoseconds\n", t0, time);
  printf("\n\nMaximum systematic error in amplitude \t= %lf %%", aSystError);
  printf("\nMaximum systematic error for timing \t= %lf nanoseconds\n\n", tSystError);
  return 0;
}

void setFileName(char *fName, int start, int N, double tau, double fs)
{
  sprintf(fName, "PFVectors/start%dN%dtau%.ffs%.f.txt", start, N, tau, fs);
  //  printf("\nfilename: %s\n", fName);
}

