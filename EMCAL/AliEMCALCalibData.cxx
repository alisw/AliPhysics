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

//_________________________________________________________________________
///*-- Author: Yves Schutz (SUBATECH)
//           : Aleksei Pavlinov (WSU); Jun 30, 2006 - ALICE numbering scheme
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for EMCAL calibration                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>

#include "AliEMCALCalibData.h"

ClassImp(AliEMCALCalibData)

//________________________________________________________________
AliEMCALCalibData::AliEMCALCalibData()
{
  // Default constructor
  Reset();
}

//________________________________________________________________
AliEMCALCalibData::AliEMCALCalibData(const char* name)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliEMCALCalibData::AliEMCALCalibData(const AliEMCALCalibData& calibda) :
  TNamed(calibda)
{
  // copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();

  Int_t nSMod = 12;
  Int_t nCol  = 48;
  Int_t nRow  = 24;
  Int_t nRow2 = 12; //Modules 11 and 12 are half modules

  for(Int_t supermodule=0; supermodule<nSMod; supermodule++) {
    if(supermodule > 10)
      nRow = nRow2;
    for(Int_t column=0; column<nCol; column++) {
      for(Int_t row=0; row<nRow; row++) {
	fADCchannel[supermodule][column][row] = 
	  calibda.GetADCchannel(supermodule,column,row);
	fADCpedestal[supermodule][column][row] = 
	  calibda.GetADCpedestal(supermodule,column,row);
      }
    }
  }
}


//________________________________________________________________
AliEMCALCalibData &AliEMCALCalibData::operator =(const AliEMCALCalibData& calibda)
{
  // assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();

  Int_t nSMod = 12;
  Int_t nCol  = 48;
  Int_t nRow  = 24;
  Int_t nRow2 = 12; //Modules 11 and 12 are half modules

  for(Int_t supermodule=0; supermodule<nSMod; supermodule++) {
    if(supermodule > 10)
      nRow = nRow2;
    for(Int_t column=0; column<nCol; column++) {
      for(Int_t row=0; row<nRow; row++) {
	fADCchannel[supermodule][column][row] = 
	  calibda.GetADCchannel(supermodule,column,row);
	fADCpedestal[supermodule][column][row] = 
	  calibda.GetADCpedestal(supermodule,column,row);
      }
    }
  }
  return *this;
}

//________________________________________________________________
AliEMCALCalibData::~AliEMCALCalibData()
{
  // Destructor
}

//________________________________________________________________
void AliEMCALCalibData::Reset()
{
  // Set all pedestals to 0 and all ADC channels widths to 1
  memset(fADCchannel ,1,12*48*24*sizeof(Float_t));
  memset(fADCpedestal,0,12*48*24*sizeof(Float_t));
}

//________________________________________________________________
void  AliEMCALCalibData::Print(Option_t *option) const
{
  // Print tables of pedestals and ADC channels widths

  Int_t nSMod = 12;
  Int_t nCol  = 48;
  Int_t nRow  = 24;
  Int_t nRow2 = 12; //Modules 11 and 12 are half modules
  if (strstr(option,"ped")) {
    printf("\n	----	Pedestal values	----\n\n");
    for (Int_t supermodule=0; supermodule<nSMod; supermodule++){
      if(supermodule > 10)
	nRow = nRow2;
      printf("============== Supermodule %d\n",supermodule+1);
      for (Int_t column=0; column<nCol; column++){
	for (Int_t row=0; row<nRow; row++){
	  printf("%4.1f",fADCpedestal[supermodule][column][row]);
	}
	printf("\n");
      }
    } 
  }


  if (strstr(option,"gain")) {
    printf("\n	----	ADC channel values	----\n\n");
    for (Int_t supermodule=0; supermodule<nSMod; supermodule++){
      if(supermodule > 10) 
	nRow = nRow2;
      printf("============== Supermodule %d\n",supermodule+1);
      for (Int_t column=0; column<nCol; column++){
	for (Int_t row=0; row<nRow; row++){
	  printf("%4.1f",fADCchannel[supermodule][column][row]);
	}
	printf("\n");
      }
    }   
  }
}

//________________________________________________________________
Float_t AliEMCALCalibData::GetADCchannel(Int_t supermodule, Int_t column, Int_t row) const
{
  // Set ADC channel witdth values
  // All indexes start from 0!
  // Supermodule, column,raw should follow the ALICE convention:
  // supermodule 0:11, column 0:47, row 0:23

  return fADCchannel[supermodule][column][row];
}

//________________________________________________________________
Float_t AliEMCALCalibData::GetADCpedestal(Int_t supermodule, Int_t column, Int_t row) const
{
  // Get ADC pedestal values
 return fADCpedestal[supermodule][column][row];
}

//________________________________________________________________
void AliEMCALCalibData::SetADCchannel(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{ 
  // Set ADC channel width values
  fADCchannel[supermodule][column][row] = value;
}

//________________________________________________________________
void AliEMCALCalibData::SetADCpedestal(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{
  // Set ADC pedestal values
  fADCpedestal[supermodule][column][row] = value;
}
