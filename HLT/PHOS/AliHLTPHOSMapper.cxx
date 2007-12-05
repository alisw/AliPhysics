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


AliHLTPHOSMapper::AliHLTPHOSMapper() : AliHLTPHOSBase(), fHw2geomapPtr(0)
{
  //  printf("\nCreating new mapper\n");
  InitAltroMapping(); 
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

  char filename[256];
  char *base =  getenv("ALICE_ROOT");
  int nChannels = 0;
  int maxaddr = 0;
  int tmpHwaddr = 0;
  int tmpZRow = 0;
  int tmpXCol = 0;
  int tmpGain = 0;

  if(base !=0)
    {
      sprintf(filename,"%s/PHOS/mapping/RCU0.data", base);

      FILE *fp = fopen(filename, "r");
      if(fp != 0)
	{
	  cout << "mapping file found" << endl;
	  fscanf(fp, "%d", &nChannels);
	  fscanf(fp, "%d", &maxaddr);
	  printf("nChannels = %d", nChannels);
	  printf("maxaddr = %d", maxaddr);
	  fHw2geomapPtr = new fAltromap[maxaddr +1]; 

	  for(int i=0; i< maxaddr + 1 ; i ++)
	    {
	      fHw2geomapPtr[i].fXCol = 0;
	      fHw2geomapPtr[i].fZRow = 0;
	      fHw2geomapPtr[i].fGain = 0;
	    }

	  for(int i=0; i<nChannels; i ++)
	    {
	      fscanf(fp, "%d %d %d %d\n", &tmpHwaddr, &tmpXCol, &tmpZRow,  &tmpGain);
	      fHw2geomapPtr[tmpHwaddr].fXCol   = tmpXCol;
	      fHw2geomapPtr[tmpHwaddr].fZRow   = tmpZRow;
	      fHw2geomapPtr[tmpHwaddr].fGain  = tmpGain;
	    }
	  
	}
      else
	{
	  cout << "ERROR could not find mapping file" << endl;
	}

    }
  else
    {
      printf("AliHLTPHOSMapper::InitAltroMapping(), ERROR environment ALICE_ROOT is not set, cannot find mapping file");
    }

} 

