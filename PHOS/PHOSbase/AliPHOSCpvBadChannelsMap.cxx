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
// PHOS Cpv bad channels map.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include "AliPHOSCpvBadChannelsMap.h"
#include "AliPHOSGeometry.h"
 
ClassImp(AliPHOSCpvBadChannelsMap)
 
//________________________________________________________________

  AliPHOSCpvBadChannelsMap::AliPHOSCpvBadChannelsMap() : fBads(-1)
{
  Reset();
}

//________________________________________________________________

void AliPHOSCpvBadChannelsMap::Reset()
{
  //Set all channels as good.
  
  for(Int_t module=0; module<AliPHOSCpvParam::kNDDL; module++) 
    for(Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++) 
      for(Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++) 
	fBadChannelCpv[module][column][row] = kFALSE;
 
  fBads=0;

}

//________________________________________________________________

void AliPHOSCpvBadChannelsMap::Reset(Int_t module)
{
  //Set all channels in module as good.
  Int_t nBadsInModule=0;
  for(Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++) 
    for(Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++)
      if(fBadChannelCpv[module][column][row]){
	fBadChannelCpv[module][column][row] = kFALSE;
	nBadsInModule++;
      }
  fBads-=nBadsInModule;

}

//________________________________________________________________

AliPHOSCpvBadChannelsMap::AliPHOSCpvBadChannelsMap(const AliPHOSCpvBadChannelsMap &map):
  TObject(map),fBads(map.fBads)
{
  //Copy constructor.
  for(Int_t module=0; module<AliPHOSCpvParam::kNDDL; module++) 
    for(Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++) 
      for(Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++) 
	fBadChannelCpv[module][column][row] = map.fBadChannelCpv[module][column][row];
 
}

//________________________________________________________________

AliPHOSCpvBadChannelsMap& AliPHOSCpvBadChannelsMap::operator= (const AliPHOSCpvBadChannelsMap &map) 
{
  //Assignment operator.

  if(this != &map) {
    fBads = map.fBads;
    for(Int_t module=0; module<AliPHOSCpvParam::kNDDL; module++) 
      for(Int_t column=0; column<AliPHOSCpvParam::kPadPcX; column++) 
	for(Int_t row=0; row<AliPHOSCpvParam::kPadPcY; row++) 
	  fBadChannelCpv[module][column][row] = map.fBadChannelCpv[module][column][row]; 
  }

  return *this;
}

//_________________________________________________________________

void AliPHOSCpvBadChannelsMap::SetBadChannel(Int_t module, Int_t col, Int_t row)
{
  // Declare a channel (module,col,row) as a bad, if it was not set before

  if (!fBadChannelCpv[module-1][col-1][row-1]) {
    fBadChannelCpv[module-1][col-1][row-1] = kTRUE;
    ++fBads; 
  }
}
//_________________________________________________________________

void AliPHOSCpvBadChannelsMap::BadChannelIds(Int_t *badIds)
{
  //Fill array badIds by the Ids of bad channels.
  //Array badIds of length GetNumOfBadChannels() should be prepared in advance. 

  if(!badIds) return;
  if(fBads <=0 ) return;

  AliPHOSGeometry* geom = AliPHOSGeometry::GetInstance();

  if(!geom)
    geom = AliPHOSGeometry::GetInstance("IHEP");

  Int_t absId;
  Int_t relId[4];

  Int_t iBad = 0;
  relId[1] =  -1; // CPV pad

  for(Int_t mod=1; mod<=AliPHOSCpvParam::kNDDL; mod++) { 
    for(Int_t col=1; col<=AliPHOSCpvParam::kPadPcX; col++) { 
      for(Int_t row=1; row<=AliPHOSCpvParam::kPadPcY; row++) {
	if(IsBadChannel(mod,col,row)) {
	  relId[0] = mod;
	  relId[2] = col;
	  relId[3] = row;
	  geom->RelToAbsNumbering(relId,absId);
	  badIds[iBad]=absId;
	  iBad++;
	}
      }
    }
  }

}
