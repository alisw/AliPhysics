#include "AliHLTAltroData.h"



AliHLTAltroData::AliHLTAltroData(): fBunchCounter(0), fBunchData(0), fData(0), fDataSize(0), fWc(0), fHadd(0), fIsComplete(true)
{


}


AliHLTAltroData::~AliHLTAltroData()
{


}


bool
AliHLTAltroData::NextBunch(AliHLTAltroBunch &altroBunch)
{
  if(fIsComplete == true)
    {

      /*
      cout <<"AliHLTAltroData::NextBunch" << endl;

      cout <<  "fDataSize = " << fDataSize << endl;
      cout <<  "fWc = " << fWc << endl;
      */

      if(fBunchCounter == 0)
	{
	  fBunchData = &fData[fDataSize - 1];
	}
      else
	{
	  //     fBunchData ++;
	}

      if(fWc < fDataSize)
	{
	  if(*fBunchData == 0){ fWc += 1;};
	  fWc += *fBunchData;
	  altroBunch.fData = fData - *fBunchData -1; ;
	  altroBunch.fBunchSize = *fBunchData -2;
	  fBunchData --;
	  altroBunch.fEndTimeBin = *fBunchData;
	  //      fBunchData = fBunchData - 

	  //  fBunchData --;
	  //  fBunchData = fBunchData  -  (altroBunch.fBunchSize +3);
	  fBunchData = fBunchData  -  (altroBunch.fBunchSize +1);

	  fBunchCounter ++;
	  //   fBunchData --;
	  //   fBunchData --;
	  //     cout <<"AliHLTAltroData::NextBunch, TRUE" << endl;
	  return true;
	}
      else
	{
	  fBunchCounter = 0;
	  //    cout <<"AliHLTAltroData::NextBunch, FALSE" << endl;
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

void
AliHLTAltroData::Reset()
{
   fWc = 0;
   fBunchCounter = 0;
}


int
AliHLTAltroData::GetChannel()
{
 return  fHadd & 0xf;
}

int
AliHLTAltroData::GetChip()
{
 return  (fHadd & 0x70) >> 4 ;
}

int
AliHLTAltroData::GetCard()
{
 return   (fHadd & 0x780) >>  7;
}


int
AliHLTAltroData::GetBranch()
{
 return   (fHadd & 0x800 ) >> 11;
}
