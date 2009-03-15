// $Id$

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

#include "AliHLTPHOSMapper.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTLogging.h"
#include "Rtypes.h"
#include "unistd.h"
#include <iostream>

using namespace PhosHLTConst;
AliHLTPHOSMapper::AliHLTPHOSMapper():  AliHLTLogging(), 
				       fHw2geomapPtr(0),
				       fIsInitializedMapping(false),
				       fSpecificationMapPtr(0)
{
  InitAltroMapping(); 
  InitDDLSpecificationMapping();
}


AliHLTPHOSMapper::~AliHLTPHOSMapper()
{
  delete []  fHw2geomapPtr;
  fHw2geomapPtr = 0;
}

void
AliHLTPHOSMapper::InitAltroMapping()
{
  // Loads mapping between Altro addresses and geometrical addresses from file

  //  char filename[256];
  char *base =  getenv("ALICE_ROOT");
  int nChannels = 0;
  int maxaddr = 0;
  int tmpHwaddr = 0;
  int tmpZRow = 0;
  int tmpXCol = 0;
  int tmpGain = 0;
  int res = 0; //OD to get rid of warnings
  if(base !=0)
    {
      sprintf(fFilepath,"%s/PHOS/mapping/RCU0.data", base);
      
      FILE *fp = fopen(fFilepath, "r");
      if(fp != 0)
	{
	  res = fscanf(fp, "%d", &nChannels);
	  res = fscanf(fp, "%d", &maxaddr);
	  fHw2geomapPtr = new fAltromap[maxaddr +1]; 

	  for(int i=0; i< maxaddr + 1 ; i ++)
	    {
	      fHw2geomapPtr[i].fXCol = 0;
	      fHw2geomapPtr[i].fZRow = 0;
	      fHw2geomapPtr[i].fGain = 0;
	    }

	  for(int i=0; i<nChannels; i ++)
	    {
	      res = fscanf(fp, "%d %d %d %d\n", &tmpHwaddr, &tmpXCol, &tmpZRow,  &tmpGain);
	      fHw2geomapPtr[tmpHwaddr].fXCol   = tmpXCol;
	      fHw2geomapPtr[tmpHwaddr].fZRow   = tmpZRow;
	      fHw2geomapPtr[tmpHwaddr].fGain  = tmpGain;
	    }
	  fIsInitializedMapping = true;	  
	  fclose(fp);
	}
      else
	{
	  fIsInitializedMapping = false;	  
	}
    }
  else
    {
      fIsInitializedMapping = false;
    }
} 

void
AliHLTPHOSMapper::InitDDLSpecificationMapping()
{
  fSpecificationMapPtr = new fDDLSpecificationMap[PhosHLTConst::NMODULES*PhosHLTConst::NRCUSPERMODULE];
  for(Int_t ddl = 0; ddl < PhosHLTConst::NMODULES*PhosHLTConst::NRCUSPERMODULE; ddl++)
    {
      
      fSpecificationMapPtr[ddl].fModId = ddl/PhosHLTConst::NRCUSPERMODULE;
      
      if(ddl%4 == 0)
	{
	  fSpecificationMapPtr[ddl].fRcuX = 0; 
	  fSpecificationMapPtr[ddl].fRcuZ = 0;
	}
      
      else if(ddl%4 == 1)
	{
	  fSpecificationMapPtr[ddl].fRcuX = 0; 
	  fSpecificationMapPtr[ddl].fRcuZ = 1;
	}
      
      else if( ddl%4 == 2)
	{
	  fSpecificationMapPtr[ddl].fRcuX = 1; 
	  fSpecificationMapPtr[ddl].fRcuZ = 0;
	}
      else
	{
	  fSpecificationMapPtr[ddl].fRcuX = 1; 
	  fSpecificationMapPtr[ddl].fRcuZ = 1;
	}
      
      fSpecificationMapPtr[ddl].fRcuZOffset = NZROWSRCU*(fSpecificationMapPtr[ddl].fRcuZ);
      fSpecificationMapPtr[ddl].fRcuXOffset = NXCOLUMNSRCU*(fSpecificationMapPtr[ddl].fRcuX);
    }
}


bool 
AliHLTPHOSMapper::GetIsInitializedMapping()
{
  return  fIsInitializedMapping;
}


char* 
AliHLTPHOSMapper::GetFilePath()
{
  return  fFilepath;
}

UShort_t
AliHLTPHOSMapper::GetChannelID(Int_t specification, Int_t hwAddress)
{
  
  Short_t index = 0;

  if(specification == 0x00001) index = 0;
  else if(specification == 0x00002) index = 1;
  else if(specification == 0x00004) index = 2;
  else if(specification == 0x00008) index = 3;

  else if(specification == 0x00010) index = 4;
  else if(specification == 0x00020) index = 5;
  else if(specification == 0x00040) index = 6;
  else if(specification == 0x00080) index = 7;

  else if(specification == 0x00100) index = 8;
  else if(specification == 0x00200) index = 9;
  else if(specification == 0x00400) index = 10;
  else if(specification == 0x00800) index = 11;

  else if(specification == 0x01000) index = 12;
  else if(specification == 0x02000) index = 13;
  else if(specification == 0x04000) index = 14;
  else if(specification == 0x08000) index = 15;

  else if(specification == 0x10000) index = 16;
  else if(specification == 0x20000) index = 17;
  else if(specification == 0x40000) index = 18;
  else if(specification == 0x80000) index = 19;

  else HLTError("Specification 0x%X not consistent with single DDL in PHOS", specification);

  return ((fHw2geomapPtr[hwAddress].fXCol + fSpecificationMapPtr[index].fRcuXOffset) |
	  ((fHw2geomapPtr[hwAddress].fZRow + fSpecificationMapPtr[index].fRcuZOffset) << 6) |
	  (fHw2geomapPtr[hwAddress].fGain << 12) |
	  fSpecificationMapPtr[index].fModId << 13);
}

void
AliHLTPHOSMapper::GetChannelCoord(UShort_t channelId, UShort_t* channelCoord)
{
  channelCoord[0] = channelId&0x3f;
  channelCoord[1] = (channelId >> 6)&0x3f;
  channelCoord[2] = (channelId >> 12)&0x1;
  channelCoord[3] = (channelId >> 13)&0x1f;
}
