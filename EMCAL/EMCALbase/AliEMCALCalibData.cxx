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

#include "AliEMCALCalibData.h"

/// \cond CLASSIMP
ClassImp(AliEMCALCalibData) ;
/// \endcond

///
/// Default constructor
//____________________________________
AliEMCALCalibData::AliEMCALCalibData() : 
TNamed(), fADCchannelRef(0)
{
  Reset();
}

///
/// Constructor
//____________________________________________________
AliEMCALCalibData::AliEMCALCalibData(const char* name) : 
TNamed(name,name),fADCchannelRef(0)
{
  Reset();
}

///
/// Copy constructor
//____________________________________________________________________
AliEMCALCalibData::AliEMCALCalibData(const AliEMCALCalibData& calibda) :
TNamed(calibda), fADCchannelRef(calibda.fADCchannelRef)
{
  SetName (calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules;
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;    //48
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;    //24
  
  for(Int_t supermodule = 0; supermodule < nSMod; supermodule++)
  {    
    // Init all SM equally, even the channels known to not exist.
    
    for(Int_t column = 0; column<nCol; column++)
    {
      for(Int_t row = 0; row<nRow; row++)
      {
        SetADCchannel      (supermodule,column,row, calibda.GetADCchannel      (supermodule,column,row));
        SetADCchannelOnline(supermodule,column,row, calibda.GetADCchannelOnline(supermodule,column,row));
        SetADCchannelDecal (supermodule,column,row, calibda.GetADCchannelDecal (supermodule,column,row));
        SetADCpedestal     (supermodule,column,row, calibda.GetADCpedestal     (supermodule,column,row));
        SetTimeChannelDecal(supermodule,column,row, calibda.GetTimeChannelDecal(supermodule,column,row));
        for(Int_t bc = 0; bc < 4; bc++)
          SetTimeChannel(supermodule,column,row,bc,calibda.GetTimeChannel(supermodule,column,row,bc));
      } // col
    } // row
  } // SM
}

///
/// Assignment operator
//________________________________________________________________________________
AliEMCALCalibData &AliEMCALCalibData::operator =(const AliEMCALCalibData& calibda)
{
  SetName (calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  
  fADCchannelRef = calibda.GetADCchannelRef() ;
  
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules;
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;    //48
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;    //24
  
  for(Int_t supermodule = 0; supermodule < nSMod; supermodule++)
  {
    // Init all SM equally, even the channels known to not exist.
    
    for(Int_t column = 0; column<nCol; column++)
    {
      for(Int_t row = 0; row<nRow; row++)
      {
        SetADCchannel      (supermodule,column,row, calibda.GetADCchannel      (supermodule,column,row));
        SetADCchannelOnline(supermodule,column,row, calibda.GetADCchannelOnline(supermodule,column,row));
        SetADCchannelDecal (supermodule,column,row, calibda.GetADCchannelDecal (supermodule,column,row));
        SetADCpedestal     (supermodule,column,row, calibda.GetADCpedestal     (supermodule,column,row));
        SetTimeChannelDecal(supermodule,column,row, calibda.GetTimeChannelDecal(supermodule,column,row));
        for(Int_t bc = 0; bc < 4; bc++)
          SetTimeChannel(supermodule,column,row,bc,calibda.GetTimeChannel(supermodule,column,row,bc));
      } // col
    } // row
  } // col
  
  return *this;
}

///
/// Set all pedestals to 0 and all ADC channels widths to 1
//_____________________________
void AliEMCALCalibData::Reset()
{  
  fADCchannelRef = 0.0162;	
  
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules; 
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;    //48
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;    //24

   for (Int_t supermodule=0; supermodule<nSMod; supermodule++)
  {
    //Init all SM equally, even the channels known to not exist.

    for (Int_t column=0; column < nCol; column++)
    {
      for (Int_t row = 0; row < nRow; row++)
      {
        SetADCchannel      (supermodule,column,row, fADCchannelRef);
        SetADCchannelOnline(supermodule,column,row, fADCchannelRef);
        SetADCchannelDecal (supermodule,column,row, 1);
        SetADCpedestal     (supermodule,column,row, 0);
        SetTimeChannelDecal(supermodule,column,row, 0);
        for(Int_t bc = 0; bc < 4; bc++)
          SetTimeChannel(supermodule,column,row, bc, 0);

      }
    }
  }	
}

///
/// Print tables of pedestals and ADC channels widths
/// options are: "gain", "ped", "online", "decal", "time", "timdecal", "all"
//____________________________________________________
void  AliEMCALCalibData::Print(Option_t *option) const
{  
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules;
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;    //48
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;    //24
  
  for (Int_t supermodule = 0; supermodule < nSMod; supermodule++)
  {
    //Init all SM equally, even the channels known to not exist.

    printf("============== Supermodule %d\n",supermodule+1);
    for (Int_t column = 0; column < nCol; column++)
    {
      for (Int_t row = 0; row < nRow; row++)
      {
        printf("[col %d,row %d] ",column, row);
        if (strstr(option,"gain") || strstr(option,"all"))
          printf("calib=%2.4f ",GetADCchannel(supermodule,column,row));
        
        if (strstr(option,"online") || strstr(option,"all"))
          printf("calib0=%2.4f ", GetADCchannelOnline(supermodule,column,row));
        
        if (strstr(option,"decal") || strstr(option,"all"))
          printf("calibDecal=%2.4f ",GetADCchannelDecal(supermodule,column,row));
        
        if (strstr(option,"ped") || strstr(option,"all"))
          printf("ped=%2.4f ", GetADCpedestal(supermodule,column,row));
        
        if (strstr(option,"time") || strstr(option,"all"))
          printf("time::bc0 =%2.4f, bc1=%2.4f, bc2=%2.4f, bc3=%2.4f ",
                 GetTimeChannel(supermodule,column,row,0), GetTimeChannel(supermodule,column,row,1), GetTimeChannel(supermodule,column,row,2), GetTimeChannel(supermodule,column,row,3));
        
        if (strstr(option,"timdecal") || strstr(option,"all"))
          printf("timeDecal=%2.4f ", GetTimeChannelDecal(supermodule,column,row));
        
         if (strstr(option,"all") || (row%4==3) ) printf("\n");
      }
    }
  }
}

///
/// Set ADC channel witdth values
/// All indexes start from 0!
/// Supermodule, column,raw should follow the ALICE convention:
/// supermodule 0:11, column 0:47, row 0:23
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//________________________________________________________________________________________
Float_t AliEMCALCalibData::GetADCchannel(Int_t supermodule, Int_t column, Int_t row) const
{  
  if(supermodule < 12) return fADCchannel    [supermodule]   [column][row];
  else                 return fADCchannelDCAL[supermodule-12][column][row];
}

///
/// Set ADC channel witdth values, first online calibration parameter
/// All indexes start from 0!
/// Supermodule, column,raw should follow the ALICE convention:
/// supermodule 0:11, column 0:47, row 0:23
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//______________________________________________________________________________________________
Float_t AliEMCALCalibData::GetADCchannelOnline(Int_t supermodule, Int_t column, Int_t row) const
{  
  if(supermodule < 12) return fADCchannelOnline    [supermodule]   [column][row];
  else                 return fADCchannelOnlineDCAL[supermodule-12][column][row];
}

///
/// Set ADC channel decalibration witdth values
/// All indexes start from 0!
/// Supermodule, column,raw should follow the ALICE convention:
/// supermodule 0:11, column 0:47, row 0:23
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//_____________________________________________________________________________________________
Float_t AliEMCALCalibData::GetADCchannelDecal(Int_t supermodule, Int_t column, Int_t row) const
{
  if(supermodule < 12) return fADCchannelDecal    [supermodule]   [column][row];
  else                 return fADCchannelDecalDCAL[supermodule-12][column][row];
}

///
/// Get ADC pedestal values
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//_________________________________________________________________________________________
Float_t AliEMCALCalibData::GetADCpedestal(Int_t supermodule, Int_t column, Int_t row) const
{  
  if(supermodule < 12) return fADCpedestal    [supermodule]   [column][row];
  else                 return fADCpedestalDCAL[supermodule-12][column][row];
}

///
/// Set ADC channel width values
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//______________________________________________________________________________________________
void AliEMCALCalibData::SetADCchannel(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{   
  if(supermodule < 12) fADCchannel    [supermodule]   [column][row] = value;
  else                 fADCchannelDCAL[supermodule-12][column][row] = value;
}

///
/// Set ADC channel online width values
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//____________________________________________________________________________________________________
void AliEMCALCalibData::SetADCchannelOnline(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{  
  if(supermodule < 12) fADCchannelOnline    [supermodule]   [column][row] = value;
  else                 fADCchannelOnlineDCAL[supermodule-12][column][row] = value;
}

///
/// Set ADC channel width values, decalibration
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//___________________________________________________________________________________________________
void AliEMCALCalibData::SetADCchannelDecal(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{ 
  if(supermodule < 12) fADCchannelDecal    [supermodule]   [column][row] = value;
  else                 fADCchannelDecalDCAL[supermodule-12][column][row] = value;
}

///
/// Set ADC pedestal values
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//_______________________________________________________________________________________________
void AliEMCALCalibData::SetADCpedestal(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{  
  if(supermodule < 12) fADCpedestal    [supermodule]   [column][row] = value;
  else                 fADCpedestalDCAL[supermodule-12][column][row] = value;
}

//*** Do not use the following methods use instead AliEMCALCalibTime:

///
/// Set channel time values
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//___________________________________________________________________________________________________
Float_t AliEMCALCalibData::GetTimeChannel(Int_t supermodule, Int_t column, Int_t row, Int_t bc) const
{  
  if(supermodule < 12) return fTimeChannel    [supermodule]   [column][row][bc];
  else                 return fTimeChannelDCAL[supermodule-12][column][row][bc];
}

///
/// Set channel time decalibrated values
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//______________________________________________________________________________________________
Float_t AliEMCALCalibData::GetTimeChannelDecal(Int_t supermodule, Int_t column, Int_t row) const
{
  if(supermodule < 12) return fTimeChannelDecal    [supermodule]   [column][row];
  else                 return fTimeChannelDecalDCAL[supermodule-12][column][row];
}

///
/// Set time per channel and bunch crossing values
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//_________________________________________________________________________________________________________
void AliEMCALCalibData::SetTimeChannel(Int_t supermodule, Int_t column, Int_t row, Int_t bc, Float_t value)
{  
  if(supermodule < 12) fTimeChannel    [supermodule]   [column][row][bc] = value;
  else                 fTimeChannelDCAL[supermodule-12][column][row][bc] = value;
}

///
/// Set decalibrated time per channel values
/// DCal has its own array, starting from SM 12, index 0, same col-row size
//____________________________________________________________________________________________________
void AliEMCALCalibData::SetTimeChannelDecal(Int_t supermodule, Int_t column, Int_t row, Float_t value)
{
  if(supermodule < 12) fTimeChannelDecal    [supermodule]   [column][row] = value;
  else                 fTimeChannelDecalDCAL[supermodule-12][column][row] = value;
}



