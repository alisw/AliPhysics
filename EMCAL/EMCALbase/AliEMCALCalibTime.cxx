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

#include "AliEMCALCalibTime.h"

/// \cond CLASSIMP
ClassImp(AliEMCALCalibTime) ;
/// \endcond

///
/// Default constructor
//__________________________________________________________
AliEMCALCalibTime::AliEMCALCalibTime() : TNamed()
{
  Reset();
}

///
/// Constructor
//____________________________________________________
AliEMCALCalibTime::AliEMCALCalibTime(const char* name) : TNamed(name,name)
{
  Reset();
}

///
/// Copy constructor
//____________________________________________________________________
AliEMCALCalibTime::AliEMCALCalibTime(const AliEMCALCalibTime& calibti) : TNamed(calibti)
{
  SetName (calibti.GetName());
  SetTitle(calibti.GetName());
  Reset();
    
  for(Int_t supermodule = 0; supermodule < fgkECALDCALModules; supermodule++)
  {
    // Init all SM equally, even the channels known to not exist.
    
    for (Int_t column = 0; column < AliEMCALGeoParams::fgkEMCALCols; column++)
    {
      for (Int_t row = 0; row < AliEMCALGeoParams::fgkEMCALRows; row++)
      {
        SetTimeChannelDecal(supermodule,column,row, calibti.GetTimeChannelDecal(supermodule,column,row));
        for(Int_t bc = 0; bc < 4; bc++)
          SetTimeChannel(supermodule,column,row,bc,calibti.GetTimeChannel(supermodule,column,row,bc));
      } // row
    } // col
  } // SM
}

///
/// Assignment operator.
//________________________________________________________________________________
AliEMCALCalibTime &AliEMCALCalibTime::operator =(const AliEMCALCalibTime& calibti)
{
  SetName (calibti.GetName());
  SetTitle(calibti.GetName());
  Reset();
    
  for(Int_t supermodule = 0; supermodule < fgkECALDCALModules; supermodule++)
  {
    // Init all SM equally, even the channels known to not exist.
    
    for (Int_t column = 0; column < AliEMCALGeoParams::fgkEMCALCols; column++)
    {
      for (Int_t row = 0; row < AliEMCALGeoParams::fgkEMCALRows; row++)
      {
        SetTimeChannelDecal(supermodule,column,row, calibti.GetTimeChannelDecal(supermodule,column,row));
        for(Int_t bc = 0; bc < 4; bc++)
          SetTimeChannel(supermodule,column,row,bc,calibti.GetTimeChannel(supermodule,column,row,bc));
      } // row
    } // col
  } // SM
  
  return *this;
}

///
/// Set all channels to 0. 
//_____________________________
void AliEMCALCalibTime::Reset()
{
  for (Int_t supermodule=0; supermodule < fgkECALDCALModules; supermodule++)
  {    
    // Init all SM equally, even the channels known to not exist.
    
    for (Int_t column = 0; column < AliEMCALGeoParams::fgkEMCALCols; column++)
    {
      for (Int_t row = 0; row < AliEMCALGeoParams::fgkEMCALRows; row++)
      {       
        SetTimeChannelDecal(supermodule,column,row, 0);
        for(Int_t bc = 0; bc < 4; bc++)
          SetTimeChannel(supermodule,column,row, bc, 0);
      } // row
    } // col
  }	// SM
}

///
/// Print tables 
//____________________________________________________
void  AliEMCALCalibTime::Print(Option_t *option) const
{    
  for (Int_t supermodule = 0; supermodule < fgkECALDCALModules; supermodule++)
  {
    printf("============== Supermodule %d\n",supermodule);
    for (Int_t column = 0; column < AliEMCALGeoParams::fgkEMCALCols; column++)
    {
      for (Int_t row = 0; row < AliEMCALGeoParams::fgkEMCALRows; row++)
      {
        
          printf("col %d, row %d; time::bc0 =%2.4f, bc1=%2.4f, bc2=%2.4f, bc3=%2.4f, decal %2.4f\n",
                 column, row,
                 GetTimeChannel     (supermodule,column,row,0), 
                 GetTimeChannel     (supermodule,column,row,1), 
                 GetTimeChannel     (supermodule,column,row,2), 
                 GetTimeChannel     (supermodule,column,row,3),
                 GetTimeChannelDecal(supermodule,column,row));
        
      } // row
    } // col
  } // SM
}



