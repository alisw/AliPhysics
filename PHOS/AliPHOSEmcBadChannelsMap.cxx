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
#include "AliPHOSGeometry.h"
 
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

//_________________________________________________________________

void AliPHOSEmcBadChannelsMap::SetBadChannel(Int_t module, Int_t col, Int_t row)
{
  // Declare a channel (module,col,row) as a bad, if it was not set before

  if (!fBadChannelEmc[module-1][col-1][row-1]) {
    fBadChannelEmc[module-1][col-1][row-1] = kTRUE;
    ++fBads; 
  }
}
//_________________________________________________________________

void AliPHOSEmcBadChannelsMap::BadChannelIds(Int_t *badIds)
{
  //Fill array badIds by the Ids of bad channels.
  //Array badIds of length GetNumOfBadChannels() should be prepared in advance. 

  if(!badIds) return;
  if(!fBads>0) return;

  AliPHOSGeometry* geom = AliPHOSGeometry::GetInstance();

  if(!geom)
    geom = AliPHOSGeometry::GetInstance("IHEP");

  Int_t absId;
  Int_t relId[4];

  Int_t iBad = 0;
  relId[1] =  0; // EMC crystal

  for(Int_t mod=1; mod<6; mod++) { 
    for(Int_t col=1; col<57; col++) { 
      for(Int_t row=1; row<65; row++) {
	if(IsBadChannel(mod,col,row)) {
	  relId[0] = mod;
	  relId[3] = col;
	  relId[2] = row;
	  geom->RelToAbsNumbering(relId,absId);
	  badIds[iBad]=absId;
	  iBad++;
	}
      }
    }
  }

}
