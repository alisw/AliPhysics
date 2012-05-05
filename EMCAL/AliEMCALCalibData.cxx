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
//           : Add decalibration and time calibration arrays: Jul 21, 2011 (GCB)
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for EMCAL calibration                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>

#include "AliEMCALCalibData.h"

ClassImp(AliEMCALCalibData)

//__________________________________________________________
AliEMCALCalibData::AliEMCALCalibData() : 
TNamed(), fADCchannelRef(0)
{
  // Default constructor
  Reset();
}

//________________________________________________________________________
AliEMCALCalibData::AliEMCALCalibData(const char* name) : 
TNamed(name,name),fADCchannelRef(0)
{
  // Constructor
  Reset();
}

//______________________________________________________________________
AliEMCALCalibData::AliEMCALCalibData(const AliEMCALCalibData& calibda) :
TNamed(calibda), fADCchannelRef(calibda.fADCchannelRef)
{
  // copy constructor
  SetName (calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules; //12
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;    //48
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;    //24
  Int_t nRow2 = AliEMCALGeoParams::fgkEMCALRows;    //12 - Modules 11 and 12 are half modules
  // in reality they are 1/3 but leave them as 1/2

  for(Int_t supermodule = 0; supermodule < nSMod; supermodule++) {
    
    if(supermodule >= 10)
      nRow = nRow2;
    
    for(Int_t column = 0; column<nCol; column++) {
      
      for(Int_t row = 0; row<nRow; row++) {
        
        fADCchannel[supermodule][column][row] = 
        calibda.GetADCchannel(supermodule,column,row);
        
        fADCchannelDecal[supermodule][column][row] = 
        calibda.GetADCchannelDecal(supermodule,column,row);
        
        fADCpedestal[supermodule][column][row] = 
        calibda.GetADCpedestal(supermodule,column,row);
        
        fTimeChannelDecal[supermodule][column][row] = 
        calibda.GetTimeChannelDecal(supermodule,column,row);
        
        for(Int_t bc = 0; bc < 4; bc++)
          fTimeChannel[supermodule][column][row][bc] = 
          calibda.GetTimeChannel(supermodule,column,row,bc);
        
      }
    }
  }
}

//________________________________________________________________
AliEMCALCalibData &AliEMCALCalibData::operator =(const AliEMCALCalibData& calibda)
{
  // assignment operator
  SetName (calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  
  fADCchannelRef = calibda.fADCchannelRef ;
  
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules; //12
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;    //48
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;    //24
  Int_t nRow2 = AliEMCALGeoParams::fgkEMCALRows/2;  //12 - Modules 11 and 12 are half modules
  // in reality they are 1/3 but leave them as 1/2

  for(Int_t supermodule = 0; supermodule < nSMod; supermodule++) {
    
    if(supermodule >= 10)
      nRow = nRow2;
    
    for(Int_t column = 0; column<nCol; column++) {
      
      for(Int_t row = 0; row<nRow; row++) {
        
        fADCchannel[supermodule][column][row] = 
        calibda.GetADCchannel(supermodule,column,row);
        
        fADCchannelDecal[supermodule][column][row] = 
        calibda.GetADCchannelDecal(supermodule,column,row);
        
        fADCpedestal[supermodule][column][row] = 
        calibda.GetADCpedestal(supermodule,column,row);
        
        fTimeChannelDecal[supermodule][column][row] = 
        calibda.GetTimeChannelDecal(supermodule,column,row);
        
        for(Int_t bc = 0; bc < 4; bc++)
          fTimeChannel[supermodule][column][row][bc] = 
          calibda.GetTimeChannel(supermodule,column,row,bc);
        
      }
    }
  }
  
  return *this;
}

//_____________________________
void AliEMCALCalibData::Reset()
{
  // Set all pedestals to 0 and all ADC channels widths to 1
  
  fADCchannelRef = 0.0162;	
  
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules; //12
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;    //48
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;    //24
  Int_t nRow2 = AliEMCALGeoParams::fgkEMCALRows/2;  //12 - Modules 11 and 12 are half modules
  // in reality they are 1/3 but leave them as 1/2

  for (Int_t supermodule=0; supermodule<nSMod; supermodule++){
    if(supermodule >= 10)
      nRow = nRow2;
    for (Int_t column=0; column < nCol; column++){
      
      for (Int_t row = 0; row < nRow; row++){
        
        fADCpedestal     [supermodule][column][row]=0.;
        
        fADCchannelDecal [supermodule][column][row]=1.;
        fADCchannel      [supermodule][column][row]=1.;
        
        fTimeChannelDecal[supermodule][column][row]=0.;
        
        for(Int_t bc = 0; bc < 4; bc++)
          fTimeChannel[supermodule][column][row][bc]=0;
        
      }
    }
  }	
}

//____________________________________________________
void  AliEMCALCalibData::Print(Option_t *option) const
{
  // Print tables of pedestals and ADC channels widths
  // options are: "gain", "ped", "decal", "time", "all"
  
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules; //12
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;    //48
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;    //24
  Int_t nRow2 = AliEMCALGeoParams::fgkEMCALRows/2;  //12 - Modules 11 and 12 are half modules
  // in reality they are 1/3 but leave them as 1/2
  
  if (strstr(option,"ped") || strstr(option,"all")) {
    printf("\n	----	Pedestal values	----\n\n");
    for (Int_t supermodule=0; supermodule<nSMod; supermodule++){
      if(supermodule >= 10)
        nRow = nRow2;
      printf("============== Supermodule %d\n",supermodule+1);
      for (Int_t column=0; column<nCol; column++){
        for (Int_t row=0; row<nRow; row++){
          printf(" %2.4f ",fADCpedestal[supermodule][column][row]);
        }
        printf("\n");
      }
    } 
  }
  
  if (strstr(option,"gain") || strstr(option,"all")) {
    printf("\n	----	ADC channel values	----\n\n");
    for (Int_t supermodule=0; supermodule<nSMod; supermodule++){
      if(supermodule >= 10) 
        nRow = nRow2;
      printf("============== Supermodule %d\n",supermodule+1);
      for (Int_t column=0; column<nCol; column++){
        for (Int_t row=0; row<nRow; row++){
          printf(" %2.4f ",fADCchannel[supermodule][column][row]);
        }
        printf("\n");
      }
    }   
  }
  
  if (strstr(option,"adcdecal") || strstr(option,"all")) {
    printf("\n	----	ADC decalibration channel values	----\n\n");
    for (Int_t supermodule=0; supermodule<nSMod; supermodule++){
      if(supermodule >= 10) 
        nRow = nRow2;
      printf("============== Supermodule %d\n",supermodule+1);
      for (Int_t column=0; column<nCol; column++){
        for (Int_t row=0; row<nRow; row++){
          printf(" %2.4f ",fADCchannelDecal[supermodule][column][row]);
        }
        printf("\n");
      }
    }   
  }
  
  if (strstr(option,"time") || strstr(option,"all")) {
    printf("\n	----	time channel values	----\n\n");
    for (Int_t supermodule=0; supermodule<nSMod; supermodule++){
      if(supermodule >= 10) 
        nRow = nRow2;
      printf("============== Supermodule %d\n",supermodule+1);
      for (Int_t column=0; column<nCol; column++){
        for (Int_t row=0; row<nRow; row++){
          for(Int_t bc = 0; bc < 4; bc++)
            printf(" %2.4f ",fTimeChannel[supermodule][column][row][bc]);
        }
        printf("\n");
      }
    }   
  }
  
  if (strstr(option,"time") || strstr(option,"all")) {
    printf("\n	----	time decalibration channel values	----\n\n");
    for (Int_t supermodule=0; supermodule<nSMod; supermodule++){
      if(supermodule >= 10) 
        nRow = nRow2;
      printf("============== Supermodule %d\n",supermodule+1);
      for (Int_t column=0; column<nCol; column++){
        for (Int_t row=0; row<nRow; row++){
          printf(" %2.4f ",fTimeChannelDecal[supermodule][column][row]);
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
Float_t AliEMCALCalibData::GetADCchannelDecal(Int_t supermodule, Int_t column, Int_t row) const
{
  // Set ADC channel decalibration witdth values
  // All indexes start from 0!
  // Supermodule, column,raw should follow the ALICE convention:
  // supermodule 0:11, column 0:47, row 0:23
  
  return fADCchannelDecal[supermodule][column][row];
}

//________________________________________________________________
Float_t AliEMCALCalibData::GetADCpedestal(Int_t supermodule, Int_t column, Int_t row) const
{
  // Get ADC pedestal values
  return fADCpedestal[supermodule][column][row];
}

//________________________________________________________________
Float_t AliEMCALCalibData::GetTimeChannel(Int_t supermodule, Int_t column, Int_t row, Int_t bc) const
{
  // Set channel time witdth values
  return fTimeChannel[supermodule][column][row][bc];
}

//________________________________________________________________
Float_t AliEMCALCalibData::GetTimeChannelDecal(Int_t supermodule, Int_t column, Int_t row) const
{
  // Set channel time witdth values
  return fTimeChannelDecal[supermodule][column][row];
}

//________________________________________________________________
void AliEMCALCalibData::SetADCchannel(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{ 
  // Set ADC channel width values
  fADCchannel[supermodule][column][row] = value;
}

//________________________________________________________________
void AliEMCALCalibData::SetADCchannelDecal(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{ 
  // Set ADC channel width values
  fADCchannelDecal[supermodule][column][row] = value;
}

//________________________________________________________________
void AliEMCALCalibData::SetADCpedestal(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{
  // Set ADC pedestal values
  fADCpedestal[supermodule][column][row] = value;
}

//________________________________________________________________
void AliEMCALCalibData::SetTimeChannel(Int_t supermodule, Int_t column, Int_t row, Int_t bc, Float_t value)
{
  // Set ADC pedestal values
  fTimeChannel[supermodule][column][row][bc] = value;
}

//________________________________________________________________
void AliEMCALCalibData::SetTimeChannelDecal(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{
  // Set ADC pedestal values
  fTimeChannelDecal[supermodule][column][row] = value;
}



