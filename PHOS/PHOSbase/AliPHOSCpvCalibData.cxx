/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
// Updated: Sergey Evdokimov on 28 Mar 2015
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
  for(Int_t module=0; module<AliPHOSCpvParam::kNDDL; module++) {
    for(Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++) {
      for(Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++) {
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
  for(Int_t module=0; module<AliPHOSCpvParam::kNDDL; module++) {
    for(Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++) {
      for(Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++) {
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
  // Set all pedestals and all ADC channels to its dummy values.
  for(Int_t module=0; module<AliPHOSCpvParam::kNDDL; module++) {
    for(Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++) {
      for(Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++) {
	fADCpedestalCpv[module][column][row] = 0;
	fADCchannelCpv[module][column][row] = 1;
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
    for (Int_t module=0; module<AliPHOSCpvParam::kNDDL; module++){
      printf("============== Module %d\n",module+1);
      for (Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++){
	for (Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++){
	  printf("%4.1f",fADCpedestalCpv[module][column][row]);
	}
	printf("\n");
      }
    }
  }

  if (strstr(option,"gain")) {
    printf("\n	----	ADC channel values	----\n\n");
    for (Int_t module=0; module<AliPHOSCpvParam::kNDDL; module++){
      printf("============== Module %d\n",module+1);
      for (Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++){
	for (Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++){
	  printf("%4.1f",fADCchannelCpv[module][column][row]);
	}
	printf("\n");
      }
    }
  }
}

//________________________________________________________________
Float_t AliPHOSCpvCalibData::GetADCchannelCpv(Int_t module, Int_t column, Int_t row) const
{
  //Return CPV calibration coefficient
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:128.

  return fADCchannelCpv[module-1][column-1][row-1];
}

//________________________________________________________________
Float_t AliPHOSCpvCalibData::GetADCpedestalCpv(Int_t module, Int_t column, Int_t row) const
{
  //Return CPV pedestal
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:128.
  return fADCpedestalCpv[module-1][column-1][row-1];
}

//________________________________________________________________
void AliPHOSCpvCalibData::SetADCchannelCpv(Int_t module, Int_t column, Int_t row, Float_t value)
{
  //Set CPV calibration coefficient
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:128.
  fADCchannelCpv[module-1][column-1][row-1] = value;
}

//________________________________________________________________
void AliPHOSCpvCalibData::SetADCpedestalCpv(Int_t module, Int_t column, Int_t row, Float_t value)
{
  //Set CPV pedestal
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:128.
  fADCpedestalCpv[module-1][column-1][row-1] = value;
}
