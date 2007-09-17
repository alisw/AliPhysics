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
				    fIsComplete(0)
{


}



AliAltroData::~AliAltroData()
{


}


Bool_t AliAltroData::NextBunch(AliAltroBunch *altroBunch)
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
	  altroBunch->SetBunchSize(*fBunchData -2);
	  fBunchData --;
	  altroBunch->SetEndTimeBin( *fBunchData );
	  //	  altroBunch->SetStartTimeBin(*fBunchData - fBunchSize);
	  fBunchData -= (altroBunch->GetBunchSize() +1);
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
      printf("\nAliAltroData::NextBunch: WARNING, dataset is not complet. 2AAA endmarker is missing ");
      printf("\nfor branch %d, card %d, chip %d, channel %d\n",  GetBranch(), GetCard(), GetChip(), GetChannel());
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
