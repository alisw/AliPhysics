/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

//_________________________________________________________________________
//
//
//   ZDC digit = ADC Channels for each PM 
//
//_________________________________________________________________________

#include "AliZDCDigit.h"


ClassImp(AliZDCDigit)

//____________________________________________________________________________
  AliZDCDigit::AliZDCDigit() 
{
  // Default constructor 
  
  fDetector = 0;
  fQuadrant = 0;
  fADCValue = 0;  
}

//____________________________________________________________________________
AliZDCDigit::AliZDCDigit(Int_t Det, Int_t Quad, Float_t ADCValue) 
{  
  // Constructor 
 
  fDetector = Det;
  fQuadrant = Quad;
  fADCValue = ADCValue;  
}

//____________________________________________________________________________
AliZDCDigit::AliZDCDigit(const AliZDCDigit & digit) 
{
  // Copy constructor

  fDetector = digit.fDetector;           
  fQuadrant = digit.fQuadrant;           
  fADCValue = digit.fADCValue;             

}
