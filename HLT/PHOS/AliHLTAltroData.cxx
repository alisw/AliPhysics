/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author:  Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTAltroData.h"


AliHLTAltroData::AliHLTAltroData(): fData(0),
				    fBunchData(0),
				    fDataSize(0),
				    fWc(0),
				    fHadd(0),
				    fBunchCounter(0),
				    fIsComplete(0)
{
  //comment

}



AliHLTAltroData::~AliHLTAltroData()
{
  //comment

}



bool
AliHLTAltroData::NextBunch(AliHLTAltroBunch *altroBunch)
{
  //comment
  if(fIsComplete == true)
    {

      if(fBunchCounter == 0)
	{
	  fBunchData = &fData[fDataSize - 1];
	  altroBunch->fData = &fData[fDataSize - 1];
	}

      if(fWc < fDataSize)
	{
	  fWc += *fBunchData;
	  altroBunch->fBunchSize = *fBunchData;
	  altroBunch->fBunchDataSize = altroBunch->fBunchSize  -2;

	  fBunchData --;
	  altroBunch->fEndTimeBin = *fBunchData;
	  fBunchData ++;

	  fBunchData = fBunchData  -  (altroBunch->fBunchSize);
	  altroBunch->fData = altroBunch->fData -  (altroBunch->fBunchSize);
	  //	  altroBunch->fData = fBunchData -  (altroBunch->fBunchSize);
	  
	  //	  fBunchData ++;
	  //	  altroBunch->fData ++; 
	  


	  fBunchCounter ++;
	  return true;

	}
      else
	{
	  fBunchCounter = 0;
	  fWc = 0;
	  return false;
	}
    }
  else
    {
      printf("\nAliHLTAltroData::NextBunch: WARNING, dataset is not complet. 2AAA endmarker is missing ");
      printf("\nfor branch %d, card %d, chip %d, channel %d\n",  GetBranch(), GetCard(), GetChip(), GetChannel());
      return false;
    }

}



/*
bool
AliHLTAltroData::NextBunch(AliHLTAltroBunch *altroBunch)
{
  if(fIsComplete == true)
    {

      if(fBunchCounter == 0)
	{
	  fBunchData = &fData[fDataSize - 1];
	}

      if(fWc < fDataSize)
	{
	  
	  //	  if(*fBunchData == 0)
	  //	    {
	  //	      fWc += 1;
	  //	    }
	  
	  fWc += *fBunchData;
	  altroBunch->fData = fData - *fBunchData -1; ;
	  altroBunch->fBunchDataSize = *fBunchData -2;
	  fBunchData --;
	  altroBunch->fEndTimeBin = *fBunchData;
	  cout <<  "*fBuncchData = " << *fBunchData << endl;
	  fBunchData = fBunchData  -  (altroBunch->fBunchDataSize +1);


	  fBunchCounter ++;
	  return true;
	}
      else
	{
	  fBunchCounter = 0;
	  fWc = 0;
	  return false;
	}
    }
  else
    {
      printf("\nAliHLTAltroData::NextBunch: WARNING, dataset is not complet. 2AAA endmarker is missing ");
      printf("\nfor branch %d, card %d, chip %d, channel %d\n",  GetBranch(), GetCard(), GetChip(), GetChannel());
      return false;
    }

}
*/


void
AliHLTAltroData::Reset()
{
  //comment
   fWc = 0;
   fBunchCounter = 0;
}


int
AliHLTAltroData::GetChannel()
{
  //comment
 return  fHadd & 0xf;
}

int
AliHLTAltroData::GetChip()
{
  //comment
 return  (fHadd & 0x70) >> 4 ;
}

int
AliHLTAltroData::GetCard()
{
  //comment
 return   (fHadd & 0x780) >>  7;
}


int
AliHLTAltroData::GetBranch()
{
  //comment
 return   (fHadd & 0x800 ) >> 11;
}
