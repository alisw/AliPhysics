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

#include "AliAltroData.h"
#include "AliAltroBunch.h"

ClassImp(AliAltroData)

AliAltroData::AliAltroData(): fData(0),
			      fBunchData(0),
			      fDataSize(0),
			      fWc(0),
			      fHadd(-1),
			      fPrevHadd(-1),
			      fBunchCounter(0),
			      fIsComplete(0),
			      fBufferLeft(0)
{


}



AliAltroData::~AliAltroData()
{


}


/*
bool
//AliHLTAltroData::NextBunch(AliHLTAltroBunch *altroBunch)
AliAltroData::NextBunch(AliAltroBunch *altroBunch)
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




 //Bool_t AliAltroData::NextBunch(AliAltroBunch *altroBunch)
int AliAltroData::NextBunch(AliAltroBunch *altroBunch)
{

  if(fIsComplete == kTRUE)
    {
      if(fBunchCounter == 0)
	{
	  fBunchData = &fData[fDataSize - 1];
	}

      if(fWc < fDataSize)
	{
	  if(*fBunchData == 0){ fWc += 1;};
	  fWc += *fBunchData;
	  altroBunch->SetData(fData + fDataSize - fWc);

	  int tmpsize =  *fBunchData -2;  

	  //	  altroBunch->SetBunchSize(         *fBunchData -2       );
	  altroBunch->SetBunchSize( tmpsize );
	  

	  fBufferLeft -= *fBunchData;
	  //	  printf("%s:%d, bufferleft = %d \n", __FILE__,  __LINE__ , fBufferLeft);
	  fBunchData --;
	  altroBunch->SetEndTimeBin( *fBunchData );
	  
	  // Matthias Oct 10 2008: those checks are certainly a bug, first the 
	  // bunch size is subtracted from fBufferLeft ... and than once again
	  // I can understand that it should not be negative but the check as
	  // committed in revision 29090 is wrong.
	  // Effectively, this is always skipping the last bunch of the last
	  // channel.
	  //	  if( (fBufferLeft <=  7 ) || ( fBufferLeft - tmpsize)  <= 7)
	  //if( fBufferLeft - tmpsize  <= 7)
	  if( fBufferLeft < 0)
	    {
	      //	      printf("%s:%d, ERROR, attempt too access buffer outside allowed range\n",  __FILE__ ,  __LINE__ );
	      return kFALSE;
	    }
	  

	  if(fBunchCounter >0)
	    {
	      int tmpret = altroBunch->CheckConsistency();
	      
	      if(tmpret != kTRUE)
		{
		  return tmpret;
		}

	      /*
		if( altroBunch->CheckConsistency() == kFALSE)
		{
		  return kFALSE;
		}
	      */
	    }

	  //	  altroBunch->SetStartTimeBin(*fBunchData - fBunchSize);
	  fBunchData -= (altroBunch->GetBunchSize() +1);

	  // PATCH from Per Thomas Hille 250408 mke sure tha
	  // Data is consistent by cheking the start timebin, should never be negative
	  if( (int)altroBunch->GetStartTimeBin( ) < 0)
	    {
	      //	      printf("ERROR altroBunch->GetStartTimeBin( ) is  %d", (int)altroBunch->GetStartTimeBin( ) );
	      return kFALSE;
	    }

	  fBunchCounter ++;
	  return kTRUE;
	}
      else
	{
	  fBunchCounter = 0;
	  fWc = 0;
	  return kFALSE;
	}
    }
  else
    {
      //     printf("\nAliAltroData::NextBunch: WARNING, dataset is not complet. 2AAA endmarker is missing ");
      //     printf("\nfor branch %d, card %d, chip %d, channel %d\n",  GetBranch(), GetCard(), GetChip(), GetChannel());
      return kFALSE;
    }

}



void AliAltroData::Reset()
{
   fWc = 0;
   fBunchCounter = 0;
}


Int_t AliAltroData::GetChannel() const
{
 return  fHadd & 0xf;
}

Int_t AliAltroData::GetChip() const
{
 return  (fHadd & 0x70) >> 4 ;
}

Int_t AliAltroData::GetCard() const
{
 return   (fHadd & 0x780) >>  7;
}


Int_t AliAltroData::GetBranch() const
{
 return   (fHadd & 0x800 ) >> 11;
}
