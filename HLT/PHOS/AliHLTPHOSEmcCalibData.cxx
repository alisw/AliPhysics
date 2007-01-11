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
 *                                                                        *
 * 7/1- 2007                                                              *
 * Author: Per Thomas Hille                                               *
 * code modified to include the two differnt gains of each detection      *
 * channel and also to include the Pedestals loaded into the ALTRO        * 
 * registers for basline subtraction and zero supression                  * 
 *                                                                        *
 *                                                                        * 
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for PHOS EmCal calibration                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTPHOSEmcCalibData.h"

ClassImp(AliHLTPHOSEmcCalibData)

//________________________________________________________________
  AliHLTPHOSEmcCalibData::AliHLTPHOSEmcCalibData(): TNamed()
{
  // Default constructor
  Reset();
}

//________________________________________________________________
AliHLTPHOSEmcCalibData::AliHLTPHOSEmcCalibData(const char* name)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliHLTPHOSEmcCalibData::AliHLTPHOSEmcCalibData(const AliHLTPHOSEmcCalibData& calibda) :
  TNamed(calibda)
{
  // copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());


  //  for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){
  for(Int_t module=0; module<PHOS_MODULES; module++) {
    for(Int_t column=0; column< PHOS_COLUMNS; column++) {
      for(Int_t row=0; row<PHOS_ROWS; row++) {
	for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){
	  fADCchannelEnergy[module][column][row][gain] =  calibda.fADCchannelEnergy[module][column][row][gain];
	  fADCpedestalEmcMeasured[module][column][row][gain] = calibda.fADCpedestalEmcMeasured[module][column][row][gain];
	}
      }  
    }
  }
}

//________________________________________________________________
AliHLTPHOSEmcCalibData &AliHLTPHOSEmcCalibData::operator =(const AliHLTPHOSEmcCalibData& calibda)
{
  // assignment operator

  if(this != &calibda) { 

    SetName(calibda.GetName());
    SetTitle(calibda.GetName());
    //    for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){
    for(Int_t module=0; module<PHOS_MODULES; module++) {
      for(Int_t column=0; column< PHOS_COLUMNS; column++) {
	for(Int_t row=0; row<PHOS_ROWS; row++) {
	  for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){
	    fADCchannelEnergy[module][column][row][gain] = calibda.fADCchannelEnergy[module][column][row][gain];
	    fADCpedestalEmcMeasured[module][column][row][gain] = calibda.fADCpedestalEmcMeasured[module][column][row][gain];
	  }
	}
      }
    }
  }
  return *this;
}

//________________________________________________________________
AliHLTPHOSEmcCalibData::~AliHLTPHOSEmcCalibData()
{
  // Destructor
}

//________________________________________________________________
void AliHLTPHOSEmcCalibData::Reset()
{
  // Set all pedestals and all ADC channels to its ideal values = 1.

  //  for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){
    for (Int_t module=0; module<PHOS_MODULES; module++){
      for (Int_t column=0; column< PHOS_COLUMNS; column++){
	for (Int_t row=0; row<PHOS_ROWS; row++){
	  for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){
	  fADCpedestalEmcMeasured[module][column][row][gain] = 0.;
	  fADCchannelEnergy[module][column][row][gain]  = 1.;
	}
      }
    }
  }
}

//________________________________________________________________
void  AliHLTPHOSEmcCalibData::Print(Option_t *option) const
{
  // Print tables of pedestals and ADC channels

  if (strstr(option,"ped")) {
    printf("\n	----	EMC Pedestal values	----\n\n");
    //    for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){
      for (Int_t module=0; module<PHOS_MODULES; module++){
	printf("============== Module %d\n",module+1);
	for (Int_t column=0; column< PHOS_COLUMNS; column++){
	  for (Int_t row=0; row<PHOS_ROWS; row++){
	    for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){
	    printf("%4.1f",fADCpedestalEmcMeasured[module][column][row][gain]);
	  }
	  printf("\n");
	}
      }
    }
  }

  if (strstr(option,"gain")) {
    printf("\n	----	EMC ADC channel values	----\n\n");
    //   for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){ 
    for (Int_t module=0; module<PHOS_MODULES; module++){
      printf("============== Module %d\n",module+1);
      for (Int_t column=0; column< PHOS_COLUMNS; column++){
	for (Int_t row=0; row<PHOS_ROWS; row++){
	  for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){ 
	    printf("%4.1f",fADCchannelEnergy[module][column][row][gain]);
	  }
	  printf("\n");
	}
      }
    }  
  }
}


void  
AliHLTPHOSEmcCalibData::MakeADCpedestalCorrectionTable()
{
  for (Int_t module=0; module<PHOS_MODULES; module++){
    printf("============== Module %d\n",module+1);
    for (Int_t column=0; column< PHOS_COLUMNS; column++){
      for (Int_t row=0; row<PHOS_ROWS; row++){
	for(Int_t gain = 0; gain < PHOS_GAINS; gain ++){ 
	  fADCpedestalCorrectionTable[module][column][row][gain] = fADCpedestalEmcMeasured[module][column][row][gain] - fADCpedestalAltroReg[module][column][row][gain];
	  //	  printf("%4.1f",fADCchannelEnergy[module][column][row][gain]);
	  //
	}
	printf("\n");
      }
    }
  }    
}


//________________________________________________________________
Float_t AliHLTPHOSEmcCalibData::GetADCchannelEnergy(Int_t module, Int_t column, Int_t row, Int_t gain) const
{
  //Return EMC calibration coefficient
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:64

  return fADCchannelEnergy[module-1][column-1][row-1][gain];
}

//________________________________________________________________
Float_t AliHLTPHOSEmcCalibData::GetADCpedestalEmcMeasured(Int_t module, Int_t column, Int_t row, Int_t gain) const
{
  //Return EMC pedestal
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:64


  return fADCpedestalEmcMeasured[module-1][column-1][row-1][gain];
}

//________________________________________________________________
void AliHLTPHOSEmcCalibData::SetADCchannelEnergy(Int_t module, Int_t column, Int_t row, Int_t gain, Float_t value)
{
  //Set EMC calibration coefficient
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:64

  fADCchannelEnergy[module-1][column-1][row-1][gain] = value;
}

//________________________________________________________________
void AliHLTPHOSEmcCalibData::SetADCpedestalEmcMeasured(Int_t module, Int_t column, Int_t row, Int_t gain, Float_t value)
{
  //Set EMC pedestal
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:64
  fADCpedestalEmcMeasured[module-1][column-1][row-1][gain] = value;
}
