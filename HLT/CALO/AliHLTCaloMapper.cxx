// $Id: AliHLTCalorimeterMapper.cxx 34622 2009-09-04 13:22:01Z odjuvsla $

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Per Thomas Hille perthi@fys.uio.no for the ALICE DCS Project.  *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
// Mapping class fro mapping
// from hardware address to geometrical address
//
//

#include "AliHLTCaloMapper.h"
#include "AliHLTLogging.h"
#include "Rtypes.h"
#include "unistd.h"
#include <iostream>
#include "AliHLTCaloCoordinate.h"

using namespace CaloHLTConst;

//ClassImp(AliHLTCaloMapper)

AliHLTCaloMapper::AliHLTCaloMapper() :  
  AliHLTLogging(), 
  fHw2geomapPtr(0),
  fIsInitializedMapping(false),
  fSpecificationMapPtr(0)
{  
  //see header file for class documentation
}


AliHLTCaloMapper::~AliHLTCaloMapper()
{
  delete []  fHw2geomapPtr;
  fHw2geomapPtr = 0;
}


//void AliHLTCaloMapper::InitAltroMapping(){
  //Virtual base class
//}

//void AliHLTCaloMapper::InitDDLSpecificationMapping() {
  //Virtual base class
//}



bool 
AliHLTCaloMapper::GetIsInitializedMapping()
{
  return  fIsInitializedMapping;
}


char* 
AliHLTCaloMapper::GetFilePath()
{
  return  fFilepath;
}


void
//AliHLTPHOSMapper::GetChannelCoord(const UShort_t channelId,    &AliHLTPHOSCoordinate channelCoord)
AliHLTCaloMapper::ChannelId2Coordinate(const UShort_t channelId,    AliHLTCaloCoordinate &channelCoord)
{
  channelCoord.fX = channelId&0x3f;
  channelCoord.fZ = (channelId >> 6)&0x3f;
  channelCoord.fGain = (channelId >> 12)&0x1;
  channelCoord.fModuleId  = (channelId >> 13)&0x1f;
  //  printf("Channel ID: 0x%X Coordinates: x = %d, z = %d, gain = %d\n", channelId, channelCoord[0], channelCoord[1], channelCoord[2]);
}
