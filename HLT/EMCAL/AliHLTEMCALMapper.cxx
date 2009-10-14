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

#include "AliHLTEMCALConstants.h"



using namespace EmcalHLTConst;

//AliHLTCaloMapper


AliHLTEMCALMapper::AliHLTEMCALMapper(const unsigned long specification ) : AliHLTCaloMapper(specification) 
{
  InitAltroMapping(specification);
  InitDDLSpecificationMapping();
  fIsInitializedMapping = true; //CRAP PTH, must check is the initilization actually went ok
}


AliHLTEMCALMapper::~AliHLTEMCALMapper()
{

}


void 
AliHLTEMCALMapper::InitAltroMapping(const unsigned long specification )
{
  char *base =  getenv("ALICE_ROOT");
  int  nChannels = 0;
  int  maxaddr = 0;
  int  tmpHwaddr = 0;
  int tmpZRow = 0;
  int tmpXCol = 0;
  int tmpGain = 0;
  int res = 0; 
  
  if(base !=0)
    {
      sprintf(fFilepath, "%s/EMCAL/mapping/%s", base,   DDL2RcuMapFileName( GetDDLFromSpec( specification ) ) ); 
      cout << __FILE__ <<":"<< __LINE__ <<"mapping filename is " <<  fFilepath << endl;
	// sprintf(fFilepath,"%s/PHOS/mapping/RCU0.data", base);
      FILE *fp = fopen(fFilepath, "r");
      if(fp != 0)
	{
	  res = fscanf(fp, "%d\n", &nChannels);
	  res = fscanf(fp, "%d\n", &maxaddr);
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
	      
	      //	      cout << __FILE__ << __LINE__ << "  tmpHwaddr  = " << tmpHwaddr << ", tmpXCol = " << (int)tmpXCol <<  ", tmpZRow = "<< (int)tmpZRow <<  ", tmpGain= "<< (int)tmpGain << endl;
	      
	      if(tmpGain < 2)
		{
		  fHw2geomapPtr[tmpHwaddr].fXCol   = (char)tmpXCol;
		  fHw2geomapPtr[tmpHwaddr].fZRow   = (char)tmpZRow;
		  fHw2geomapPtr[tmpHwaddr].fGain  =  (char)tmpGain;
		} 
	    }
	  fIsInitializedMapping = true;	  
	  fclose(fp);
	}
      else
	{
	  cout << __FUNCTION__ << ":"<<__FILE__<<":"<< __LINE__ << "ERROR, could not open mapping file %s" <<  fFilepath << endl;
      	  fIsInitializedMapping = false;	  
	}
    }
  else
    {
      fIsInitializedMapping = false;
    }
}


void 
AliHLTEMCALMapper::InitDDLSpecificationMapping()
{
  fSpecificationMapPtr = new fDDLSpecificationMap[EmcalHLTConst::NMODULES*EmcalHLTConst::NRCUSPERMODULE];
  
  for(Int_t ddl = 0; ddl < EmcalHLTConst::NMODULES*EmcalHLTConst::NRCUSPERMODULE; ddl++)
    {
      fSpecificationMapPtr[ddl].fModId = ddl/EmcalHLTConst::NRCUSPERMODULE;
      fSpecificationMapPtr[ddl].fRcuX = 0; 
      fSpecificationMapPtr[ddl].fRcuZ = ddl%2; 
      //      fSpecificationMapPtr[ddl].fRcuZOffset = NZROWSRCU*(fSpecificationMapPtr[ddl].fRcuZ);
      //      fSpecificationMapPtr[ddl].fRcuXOffset = NXCOLUMNSRCU*(fSpecificationMapPtr[ddl].fRcuX);
    }
}



//RCU1C.data


const char* 
AliHLTEMCALMapper::DDL2RcuMapFileName(const int ddlIndex) const //0=4608, 1=4607 etc...
{
  static char rname[256];
  char tmpSide;
  
  if ( ddlIndex%NRCUSPERSECTOR <2)
    {
      tmpSide  = 'A';
    }
  else
    {
      tmpSide  = 'C';
    }
  
  sprintf(rname,"RCU%d%c.data", ddlIndex/NRCUSPERSECTOR, tmpSide );
  return rname;
  // rname.fSector = ddlIndex/NRCUSPERSECTOR;
}
 

/*
unsigned long 
AliHLTEMCALMapper::GetSpecFromDDLIndex( const int ddlindex )
{
  return ( (unsigned long)1  <<  ddlindex ));
}
*/

/*
AliHLTEMCALMapper::GlobalX2ModuleId( const int globalX )
{
  return globalX/NXCOLUMNSMOD;  
}
*/

 /*
static  const int 
AliHLTEMCALMapper::GlobalZ2ModuleId( const int globalZ )
{
  return globalZ/NZROWSMOD;
  }
 */


  /*
 const int 
 AliHLTEMCALMapper::Global2ModuleId( const int globalZ,  const int globalX )
{
  int tmpModX = 
 
}
  */
