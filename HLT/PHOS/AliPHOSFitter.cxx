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

#include "AliPHOSFitter.h"
#include <iostream>

using std::cout;
using std::endl;

AliPHOSFitter::AliPHOSFitter()
{
  cout <<"You cannot invoke the Fitter without arguments"<<endl;;
  fFloatDataPtr  = 0;
  kfMCovarPtrPtr = 0;
  fPCovarPtrPtr  = 0;
}

AliPHOSFitter::AliPHOSFitter(double *dtaPtr, double fs)
{
  fFloatDataPtr = dtaPtr;  
  fSampleFrequency = fs;
} //end AliPHOSFitter    

AliPHOSFitter::~AliPHOSFitter()
{
  delete fFloatDataPtr;
  delete kfMCovarPtrPtr;
  delete fPCovarPtrPtr;
  fFloatDataPtr  = 0;
  kfMCovarPtrPtr = 0;
  fPCovarPtrPtr  = 0;
} //end AliPHOSFitter

void 
AliPHOSFitter::BaselineCorrection(double *dataPtr, int N)
{
  cout << "Baseline correction not yet implemeted" << endl;
} //end BaselineCorrection

void 
AliPHOSFitter::BaselineCorrection(double *dataPtr, double baselineValue)
{
  cout << "Baseline correction not yet implemeted" << endl;
} //end BaslineCorrection

void 
AliPHOSFitter::FitChiSquare(int start, int lenght)
{
  cout << "Chi square fit is not yet implemented" << endl;
} //end FitChiSquare

void 
AliPHOSFitter::FitChiSquare(int start, int lenght, double tGuess, double aGues)
{
  cout << "Chi square fit is not yet implemented" << endl;
} //end FitChiSquare

void 
AliPHOSFitter::FitKLevel(double kLevel, int start, int lenght)
{
  cout << "Fitting with the k-level method is not yet implemented" << endl;
} //end FitKLevel

void 
AliPHOSFitter::FitLeastMeanSquare(int start, int lenght, const double **kmCovarPtrPtr, double **pCovarPtrPtr)
{
  cout << "Optimized least square fit is not yet implemented" << endl;
} //end FitLeastMeanSquare

void 
AliPHOSFitter::FitLeastMeanSquare(int start, int lenght, const double **kmCovarPtrPtr, double **pCovarPtrPtr, double tGuess, double aGues)
{
  cout << "Optimized least square fit is not yet implemented" << endl;
} //end FitLeastMeanSquare

void 
AliPHOSFitter::FitPeakFinder(int start, int length, double *tVector, double *aVector)
{
  fDTof = 0;
  fDAmpl = 0;

  double tmpAmpl[length];;
  double tmpTime[length];

  for(int i=0; i < length; i++)
    {  
      fDAmpl += aVector[i]*fFloatDataPtr[i];    
    }
  
  for(int i=0; i < length; i++)
    {   
      tmpTime[i] = tVector[i]*fFloatDataPtr[i];
      fDTof = fDTof + tmpTime[i]; 
   }

  fDTof = fDTof/fDAmpl;
  //thats all 
} //end FitPeakFinder

int 
AliPHOSFitter::FindStartIndex(double treshold)
{
  cout << "Find Start index not yet implemented" << endl;
} //end FindStartIndex


float
AliPHOSFitter::GetTiming()
{
  double tmp;
  //  tmp = (fDTof/fDAmpl);
  return fDTof;
  //  return tmp;
} //end GetTiming


float
AliPHOSFitter::GetEnergy()
{
  return fDAmpl;
} //end GetEnergy


void 
AliPHOSFitter::SetData(double *data)
{
  cout << "Set data not yet implemented" << endl;
}

//private section

void 
AliPHOSFitter::MakeInitialGuess()
{
  cout << "Make initial guess not yet implemeted" << endl;
}


void 
AliPHOSFitter::MakeInitialGuess(int treshold)
{
 cout << "Make initial guess not yet implemeted" << endl;  
}

