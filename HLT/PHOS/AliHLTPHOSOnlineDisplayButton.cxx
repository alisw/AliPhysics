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
#include "AliHLTPHOSOnlineDisplayButton.h"
#include <iostream>
#include "AliHLTPHOSOnlineDisplay.h"

using namespace std;


AliHLTPHOSOnlineDisplayButton::AliHLTPHOSOnlineDisplayButton()
{

}


//AliHLTPHOSOnlineDisplayButton::AliHLTPHOSOnlineDisplayButton(TGGroupFrame *gfPtr,char opt, char *name)
AliHLTPHOSOnlineDisplayButton::AliHLTPHOSOnlineDisplayButton(AliHLTPHOSOnlineDisplay *onlineDisplayPtr,  TGMainFrame *gfPtr,char opt, char *name)
 : TGTextButton(gfPtr, name)
{
  //  fOnlineDisplayPtr;
  fOnlineDisplayPtr = onlineDisplayPtr;
  command = opt;
}



AliHLTPHOSOnlineDisplayButton::~AliHLTPHOSOnlineDisplayButton()
{

}



Bool_t
AliHLTPHOSOnlineDisplayButton::HandleButton(Event_t* event)
{
  AllowStayDown(kFALSE);

  if(event->fType == kButtonPress) 
    {
      AllowStayDown(kFALSE);
      
      switch(command)
	{
	case 'r': //First get configuration comment
	  cout << "AliHLTPHOSOnlineDisplayButton::HandleButton,   getting rawdata"<< endl;  
	  //	  fOnlineDisplayPtr->ShowRawData();
	  break;
	default:
	  //	  MainGui::DisplayMessage("illegal command");
	  cout << "illegal command"  << endl;
	  break;
	}//end switch

    }//end if

}//end HandleButton
