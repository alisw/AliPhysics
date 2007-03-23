/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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
// PHOS EmCal bad channels map.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include "AliPHOSEmcBadChannelsMap.h"
 
ClassImp(AliPHOSEmcBadChannelsMap)
 
//________________________________________________________________

  AliPHOSEmcBadChannelsMap::AliPHOSEmcBadChannelsMap() : fBads(-1)
{
  Reset();
}

//________________________________________________________________

void AliPHOSEmcBadChannelsMap::Reset()
{
  //Set all channels as good.
  
  for(Int_t module=0; module<5; module++) 
    for(Int_t column=0; column<56; column++) 
      for(Int_t row=0; row<64; row++) 
	fBadChannelEmc[module][column][row] = kFALSE;
 
  fBads=0;

}

//________________________________________________________________

AliPHOSEmcBadChannelsMap::AliPHOSEmcBadChannelsMap(const AliPHOSEmcBadChannelsMap &map):
  TObject(map),fBads(map.fBads)
{
  //Copy constructor.

  for(Int_t module=0; module<5; module++) 
    for(Int_t column=0; column<56; column++) 
      for(Int_t row=0; row<64; row++) 
	fBadChannelEmc[module][column][row] = map.fBadChannelEmc[module][column][row];
 
}

//________________________________________________________________

AliPHOSEmcBadChannelsMap& AliPHOSEmcBadChannelsMap::operator= (const AliPHOSEmcBadChannelsMap &map) 
{
  //Assignment operator.

  if(this != &map) {
    fBads = map.fBads;
    for(Int_t module=0; module<5; module++) 
      for(Int_t column=0; column<56; column++) 
	for(Int_t row=0; row<64; row++) 
	  fBadChannelEmc[module][column][row] = map.fBadChannelEmc[module][column][row]; 
  }

  return *this;
}
