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

//**********************************************************************
//
//   		Class for ZDC digit 
//   	      ADC Channels for each PM 
//   	   5 for hadronic ZDCs 1 for EM ZDCs
//
//**********************************************************************

#include "AliZDCDigit.h"


ClassImp(AliZDCDigit)

//____________________________________________________________________________
  AliZDCDigit::AliZDCDigit() 
{
  // Default constructor 
  
  fSector[0]   = 0;
  fSector[1]   = 0;
  fADCValue[0] = 0;  
  fADCValue[1] = 0;  
}

//____________________________________________________________________________
AliZDCDigit::AliZDCDigit(Int_t *Sector, Int_t *ADCValue) 
{  
  // Constructor 
 
  Int_t i;
  for(i=0; i<2; i++) fSector[i] = Sector[i];
  for(i=0; i<2; i++) fADCValue[i] = ADCValue[i];  
}

//____________________________________________________________________________
AliZDCDigit::AliZDCDigit(const AliZDCDigit & digit):TObject(digit)
{
  // Copy constructor

  fSector[0]   = digit.fSector[0];           
  fSector[1]   = digit.fSector[1];           
  fADCValue[0] = digit.fADCValue[0];             
  fADCValue[1] = digit.fADCValue[1];             

}

