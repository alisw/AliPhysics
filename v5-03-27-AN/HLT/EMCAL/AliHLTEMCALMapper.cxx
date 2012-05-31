// $Id$

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
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

#include "AliHLTEMCALMapper.h"
#include "assert.h"
#include "AliHLTCaloConstants.h"
#include "TString.h"

using EMCAL::NXCOLUMNSMOD;
using EMCAL::NZROWSMOD;
using EMCAL::NMODULES;
using EMCAL::NRCUSPERMODULE;
using EMCAL::NRCUSPERSECTOR;
using EMCAL::MAXHWADDR;
using EMCAL::MAXCHANNELS; 

AliHLTEMCALMapper::AliHLTEMCALMapper(const unsigned long specification ) : AliHLTCaloMapper(specification, "EMCAL")
{
  fCellSize = 6;
  InitAltroMapping(specification);
  InitDDLSpecificationMapping();
  fIsInitializedMapping = true; //CRAP PTH, must check if the initilization actually went ok
}


AliHLTEMCALMapper::~AliHLTEMCALMapper()
{

}


void 
AliHLTEMCALMapper::GetLocalCoord(const int channelId, Float_t* localCoord) const
{
  localCoord[0] = ( ( Float_t )(channelId&0x3f) - NXCOLUMNSMOD/2)  * fCellSize;
  localCoord[1] = ( (Float_t)((channelId >> 6)&0x3f) - NZROWSMOD/2) * fCellSize;
}


Bool_t 
AliHLTEMCALMapper::InitAltroMapping(const unsigned long specification )
{
  TString base(getenv("ALICE_ROOT"));
  int  nChannels = 0;
  int  maxaddr = 0; // keep as dummy for now
  int  tmpHwaddr = 0;
  int tmpZRow = 0;
  int tmpXCol = 0;
  int tmpGain = 0;
  int res = 0; 
  
  if(!base.IsNull())
    {
      snprintf(fFilepath,FILEPATHMAXLENGTH,"%s/EMCAL/mapping/%s", base.Data(),   DDL2RcuMapFileName( GetDDLFromSpec( specification ) ) ); 
      snprintf(fFilepath, FILEPATHMAXLENGTH,"%s/EMCAL/mapping/%s", base.Data(),   DDL2RcuMapFileName( GetDDLFromSpec( specification ) ) ); 

      FILE *fp = fopen(fFilepath, "r");
      if(fp != 0)
	{
	  res = fscanf(fp, "%d\n", &nChannels);
	  res = fscanf(fp, "%d\n", &maxaddr);
	  //	  fHw2geomapPtr = new fAltromap[maxaddr +1]; 
	  fHw2geomapPtr = new fAltromap[MAXHWADDR +1];
	
 	  for(int i=0; i< MAXHWADDR + 1 ; i ++) 
	    //for(int i=0; i< maxaddr + 1 ; i ++)
	    {
	      fHw2geomapPtr[i].fXCol = 0;
	      fHw2geomapPtr[i].fZRow = 0;
	      fHw2geomapPtr[i].fGain = 0;
	    }
	  // MAXCHANNELS
	  if( nChannels > 0 && nChannels <= MAXCHANNELS )
	    {
	      for(int i=0; i<nChannels; i ++)
		{
		  res = fscanf(fp, "%d %d %d %d\n", &tmpHwaddr, &tmpXCol, &tmpZRow,  &tmpGain);
		  
		  if(tmpGain < 2)
		    {
		      if( tmpHwaddr <= MAXHWADDR  &&  tmpHwaddr >= 0  )
			{
			  fHw2geomapPtr[tmpHwaddr].fXCol   = (char)tmpXCol;
			  fHw2geomapPtr[tmpHwaddr].fZRow   = (char)tmpZRow;
			  fHw2geomapPtr[tmpHwaddr].fGain  =  (char)tmpGain;
			}
		    } 
		}
	      fIsInitializedMapping = true;	  
	    }
	  else
	    {
	      fIsInitializedMapping = false;	   
	    }
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
    
    return fIsInitializedMapping;
}


void 
AliHLTEMCALMapper::InitDDLSpecificationMapping()
{
  if (fSpecificationMapPtr) delete [] fSpecificationMapPtr;
  fSpecificationMapPtr = new fDDLSpecificationMap[NMODULES*NRCUSPERMODULE];
  
  for(Int_t ddl = 0; ddl < NMODULES*NRCUSPERMODULE; ddl++)
    {
      fSpecificationMapPtr[ddl].fModId = ddl/( NRCUSPERMODULE );
    }
}

const char* 
AliHLTEMCALMapper::DDL2RcuMapFileName(const int ddlIndex) const //0=4608, 1=4607 etc...
{
  const int rnamelen=256;
  static char rname[rnamelen];
  char tmpSide;
  
  if ( ddlIndex%NRCUSPERSECTOR <2)
    {
      tmpSide  = 'A';
    }
  else
    {
      tmpSide  = 'C';
    }
  int tmprcuindex = ddlIndex%2;
  snprintf(rname, rnamelen, "RCU%d%c.data",  tmprcuindex,  tmpSide );
  //sprintf(rname,"RCU%d%c.data", ddlIndex/NRCUSPERSECTOR, tmpSide );
  return rname;
}
 
