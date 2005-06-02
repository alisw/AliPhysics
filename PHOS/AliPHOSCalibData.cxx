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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for PHOS calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliPHOSCalibData.h"

ClassImp(AliPHOSCalibData)

//________________________________________________________________
AliPHOSCalibData::AliPHOSCalibData()
{
  // Default constructor
  Reset();
}

//________________________________________________________________
AliPHOSCalibData::AliPHOSCalibData(const char* name)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliPHOSCalibData::AliPHOSCalibData(const AliPHOSCalibData& calibda) :
  TNamed(calibda)
{
  // copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<64; column++) {
      for(Int_t row=0; row<64; row++) {
	fADCchannelEmc[module][column][row] = calibda.GetADCchannelEmc(module,column,row);
	fADCpedestalEmc[module][column][row] = calibda.GetADCpedestalEmc(module,column,row);
      }
    }
  }
}

//________________________________________________________________
AliPHOSCalibData &AliPHOSCalibData::operator =(const AliPHOSCalibData& calibda)
{
  // assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<64; column++) {
      for(Int_t row=0; row<64; row++) {
	fADCchannelEmc[module][column][row] = calibda.GetADCchannelEmc(module,column,row);
	fADCpedestalEmc[module][column][row] = calibda.GetADCpedestalEmc(module,column,row);
      }
    }
  }
  return *this;
}

//________________________________________________________________
AliPHOSCalibData::~AliPHOSCalibData()
{
  // Destructor
}

//________________________________________________________________
void AliPHOSCalibData::Reset()
{
  // Set all pedestals to 0 and all ADC channels to 1
  memset(fADCchannelEmc ,0,5*64*56*sizeof(Float_t));
  memset(fADCpedestalEmc,1,5*64*56*sizeof(Float_t));
}

//________________________________________________________________
void  AliPHOSCalibData::Print(Option_t *option) const
{
  // Print tables of pedestals and ADC channels

  if (strstr(option,"ped")) {
    printf("\n	----	Pedestal values	----\n\n");
    for (Int_t module=0; module<5; module++){
      printf("============== Module %d\n",module+1);
      for (Int_t column=0; column<64; column++){
	for (Int_t row=0; row<64; row++){
	  printf("%4.1f",fADCpedestalEmc[module][column][row]);
	}
	printf("\n");
      }
    }
  }

  if (strstr(option,"gain")) {
    printf("\n	----	ADC channel values	----\n\n");
    for (Int_t module=0; module<5; module++){
      printf("============== Module %d\n",module+1);
      for (Int_t column=0; column<64; column++){
	for (Int_t row=0; row<64; row++){
	  printf("%4.1f",fADCchannelEmc[module][column][row]);
	}
	printf("\n");
      }
    }
  }
}
