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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for CPV calibration.                                                //
// Author: Boris Polichtchouk (Boris.Polichtchouk@cern.ch).                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliPHOSCpvCalibData.h"

ClassImp(AliPHOSCpvCalibData)

//________________________________________________________________
  AliPHOSCpvCalibData::AliPHOSCpvCalibData() : TNamed()
{
  // Default constructor
  Reset();
}

//________________________________________________________________
AliPHOSCpvCalibData::AliPHOSCpvCalibData(const char* name)
{
  // Constructor
  TString namst = "CalibCPV_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliPHOSCpvCalibData::AliPHOSCpvCalibData(const AliPHOSCpvCalibData& calibda) :
  TNamed(calibda)
{
  // copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<64; column++) {
      for(Int_t row=0; row<128; row++) {
	fADCchannelCpv[module][column][row] = calibda.GetADCchannelCpv(module,column,row);
	fADCpedestalCpv[module][column][row] = calibda.GetADCpedestalCpv(module,column,row);
      }
    }
  }
}

//________________________________________________________________
AliPHOSCpvCalibData &AliPHOSCpvCalibData::operator =(const AliPHOSCpvCalibData& calibda)
{
  // assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<64; column++) {
      for(Int_t row=0; row<128; row++) {
	fADCchannelCpv[module][column][row] = calibda.GetADCchannelCpv(module,column,row);
	fADCpedestalCpv[module][column][row] = calibda.GetADCpedestalCpv(module,column,row);
      }
    }
  }
  return *this;
}

//________________________________________________________________
AliPHOSCpvCalibData::~AliPHOSCpvCalibData()
{
  // Destructor
}

//________________________________________________________________
void AliPHOSCpvCalibData::Reset()
{
  // Set all pedestals and all ADC channels to its default values.

  for (Int_t module=0; module<5; module++){
    for (Int_t column=0; column<64; column++){
      for (Int_t row=0; row<128; row++){
	fADCpedestalCpv[module][column][row] = 0.012;
	fADCchannelCpv[module][column][row] = 0.0012;
      }
    }
  }

}

//________________________________________________________________
void  AliPHOSCpvCalibData::Print(Option_t *option) const
{
  // Print tables of pedestals and ADC channels

  if (strstr(option,"ped")) {
    printf("\n	----	Pedestal values	----\n\n");
    for (Int_t module=0; module<5; module++){
      printf("============== Module %d\n",module+1);
      for (Int_t column=0; column<64; column++){
	for (Int_t row=0; row<128; row++){
	  printf("%4.1f",fADCpedestalCpv[module][column][row]);
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
	for (Int_t row=0; row<128; row++){
	  printf("%4.1f",fADCchannelCpv[module][column][row]);
	}
	printf("\n");
      }
    }
  }
}

Float_t AliPHOSCpvCalibData::GetADCchannelCpv(Int_t module, Int_t column, Int_t row) const
{
  //CPV pads
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:64, row 1:128.

  return fADCchannelCpv[module-1][column-1][row-1];
}

Float_t AliPHOSCpvCalibData::GetADCpedestalCpv(Int_t module, Int_t column, Int_t row) const
{
  return fADCpedestalCpv[module-1][column-1][row-1];
}

void AliPHOSCpvCalibData::SetADCchannelCpv(Int_t module, Int_t column, Int_t row, Float_t value)
{
  fADCchannelCpv[module-1][column-1][row-1] = value;
}

void AliPHOSCpvCalibData::SetADCpedestalCpv(Int_t module, Int_t column, Int_t row, Float_t value)
{
  fADCpedestalCpv[module-1][column-1][row-1] = value;
}
