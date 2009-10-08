// $Id: AliHLTEMCALOnlineDisplayButton.cxx 29824 2008-11-10 13:43:55Z richterm $

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
#include "AliHLTEMCALOnlineDisplayButton.h"
#include <iostream>
#include "AliHLTEMCALOnlineDisplay.h"

using namespace std;


AliHLTEMCALOnlineDisplayButton::AliHLTEMCALOnlineDisplayButton()
{

}


//AliHLTEMCALOnlineDisplayButton::AliHLTEMCALOnlineDisplayButton(TGGroupFrame *gfPtr,char opt, char *name)
AliHLTEMCALOnlineDisplayButton::AliHLTEMCALOnlineDisplayButton(AliHLTEMCALOnlineDisplay *onlineDisplayPtr,  TGMainFrame *gfPtr,char opt, char *name)
 : TGTextButton(gfPtr, name)
{
  //  fOnlineDisplayPtr;
  fOnlineDisplayPtr = onlineDisplayPtr;
  command = opt;
}



AliHLTEMCALOnlineDisplayButton::~AliHLTEMCALOnlineDisplayButton()
{

}



Bool_t
AliHLTEMCALOnlineDisplayButton::HandleButton(Event_t* event)
{
  AllowStayDown(kFALSE);

  if(event->fType == kButtonPress) 
    {
      AllowStayDown(kFALSE);
      
      switch(command)
	{
	case 'r': //First get configuration comment
	  cout << "AliHLTEMCALOnlineDisplayButton::HandleButton,   getting rawdata"<< endl;  
	  //	  fOnlineDisplayPtr->ShowRawData();
	  break;
	default:
	  //	  MainGui::DisplayMessage("illegal command");
	  cout << "illegal command"  << endl;
	  break;
	}//end switch

    }//end if

}//end HandleButton
