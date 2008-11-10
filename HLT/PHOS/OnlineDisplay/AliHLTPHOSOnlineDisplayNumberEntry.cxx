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


#include "AliHLTPHOSOnlineDisplayNumberEntry.h"
//#include "MainGui.h"

AliHLTPHOSOnlineDisplayNumberEntry::AliHLTPHOSOnlineDisplayNumberEntry()
{

}

AliHLTPHOSOnlineDisplayNumberEntry::~AliHLTPHOSOnlineDisplayNumberEntry()
{

}

AliHLTPHOSOnlineDisplayNumberEntry::AliHLTPHOSOnlineDisplayNumberEntry(const TGWindow* parent, Double_t val, Int_t digitwidth, Int_t id, 
				 TGNumberFormat::EStyle style, TGNumberFormat::EAttribute attr, 
				 TGNumberFormat::ELimit limits,
				 Double_t min, Double_t max):  
  TGNumberEntry(parent , val, digitwidth, id, style, attr, limits, min, max)
{
  //  buttonType = buttType; //c = config Id, r = readout region entry
  fButtonToNum = kFALSE;
}

void
AliHLTPHOSOnlineDisplayNumberEntry::ValueChanged(Long_t t)
{
  int tmp = GetIntNumber();
 //  printf("\nnumberentry:ValueChanged: walue gotten was: %d\n", tmp);

  if(t == 10000)
    {
      if(tmp > lowLimit)
	{
	  SetIntNumber(tmp -1);
	} 
      else
	{

	}
    }
  else if(t == 0)
    {
      if(tmp < (highLimit))
	{
	  SetIntNumber(tmp +1); 

	}
      else 
	{ 

	}
    }
}

void
AliHLTPHOSOnlineDisplayNumberEntry::ValueSet(Long_t t)
{
  int tmp = GetIntNumber();

  // printf("\nnumberentry:ValueSet: walue gotten was: %d\n", tmp);

  if(t == 10000)
    {
      if(tmp > highLimit)
	{
	  SetIntNumber(highLimit);
	  tmp = highLimit;
	}
      else if(tmp < lowLimit)
	{         
	  SetIntNumber(lowLimit); 
	  tmp = lowLimit;
	}
    }
  else if(t == 0)
    {
      if(tmp > highLimit)
	{
	  SetIntNumber(highLimit);
	  tmp = highLimit;
	}
      else if(tmp < lowLimit)
	{
	  SetIntNumber(lowLimit); 
	  
	  tmp = lowLimit;
	}
    }

  //  if(buttonType == 'c')
  //   {
  //     MainGui::GetConfigInfo(tmp); 
  //   }
}

void
AliHLTPHOSOnlineDisplayNumberEntry::SetButtonType(char c)
{
  buttonType = c;
}




void
AliHLTPHOSOnlineDisplayNumberEntry::SetLimits(int low, int high)
{
  lowLimit  = low;
  highLimit = high;
}
