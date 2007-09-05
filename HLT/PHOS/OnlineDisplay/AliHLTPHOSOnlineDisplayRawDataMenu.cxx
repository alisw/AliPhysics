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
#include "AliHLTPHOSOnlineDisplayRawDataMenu.h"
#include "TGFrame.h"
#include "TGLabel.h"
#include "AliHLTPHOSOnlineDisplayNumberEntry.h"
#include "AliHLTPHOSOnlineDisplayButton.h"
//#include "TGMainFrame.h"

//#include "TGWindow.h"

#include <iostream>

using namespace std;


AliHLTPHOSOnlineDisplayRawDataMenu::AliHLTPHOSOnlineDisplayRawDataMenu()
{
  cout << "ERROR, you cannot invoke the Online display without arguments  " << endl;
}


AliHLTPHOSOnlineDisplayRawDataMenu::AliHLTPHOSOnlineDisplayRawDataMenu(AliHLTPHOSOnlineDisplay *onlineDisplayPtr)
{
  fOnlineDisplayPtr =  onlineDisplayPtr;
  //  fEventPtr = eventPtr;
  cout << "creating new AliHLTPHOSOnlineDisplayRawDataMen "  <<endl;
  //  fWindowPtr  =  new   TGWindow(); 
  fWindowPtr  =  new   TGMainFrame(); 
  //  fFramePtr   =  new   TGGroupFrame(fWindowPtr, "HELLO WORLD") ;

  startLabelPtr = new TGLabel(fWindowPtr, "From");
  startLabelPtr-> MoveResize( 30, 140, 50, 20);
  endLabelPtr = new TGLabel(fWindowPtr, "To");
  endLabelPtr-> MoveResize(100, 140, 50, 20);
  startZInputPtr = new AliHLTPHOSOnlineDisplayNumberEntry(fWindowPtr, 0,  5, -1, (TGNumberFormat::EStyle) 5);
  //  startZInputPtr = new AliHLTPHOSOnlineDisplayNumberEntry(fWindowPtr, startZ,  5, -1, (TGNumberFormat::EStyle) 5);

  startZInputPtr->MoveResize( 30, 160, 50, 18);
  startZInputPtr->SetLimits(0, 55);
  endZInputPtr = new AliHLTPHOSOnlineDisplayNumberEntry(fWindowPtr, 0, 5, -1, (TGNumberFormat::EStyle) 5);
  endZInputPtr->MoveResize( 100,  160, 50, 18); 
  endZInputPtr->SetLimits(0, 55);
  startXInputPtr = new AliHLTPHOSOnlineDisplayNumberEntry(fWindowPtr, 0, 5, -1, (TGNumberFormat::EStyle) 5);
  startXInputPtr->MoveResize(30, 180, 50, 18);
  startXInputPtr->SetLimits(0, 63*5);
  endXInputPtr = new AliHLTPHOSOnlineDisplayNumberEntry(fWindowPtr, 0, 5, -1, (TGNumberFormat::EStyle) 5);
  endXInputPtr->MoveResize(  100, 180, 50, 18); 
  endXInputPtr->SetLimits(0, 63*5);

  gainInputPtr = new AliHLTPHOSOnlineDisplayNumberEntry(fWindowPtr, 0, 5, -1, (TGNumberFormat::EStyle) 5);
  gainInputPtr->MoveResize(  190, 180, 50, 18); 
  gainInputPtr->SetLimits(0, 1);

  gainLabelPtr = new TGLabel(fWindowPtr, "Gain");
  gainLabelPtr-> MoveResize(  190,  150, 30, 30);

  zLabelPtr = new TGLabel(fWindowPtr, "Z");
  zLabelPtr-> MoveResize(     10,   160, 20, 20);
  xLabelPtr = new TGLabel(fWindowPtr, "X");
  xLabelPtr-> MoveResize(     10,  180, 20, 20);

  fGetDataButtonPtr = new  AliHLTPHOSOnlineDisplayButton(fOnlineDisplayPtr, fWindowPtr, 'r',  "show rawdata");
  fGetDataButtonPtr->MoveResize(70,  90,  150, 20); 
  //  applyFeeButtPtr  = new PhosMenuButton(applyApdMenuPtr, 'c', "Apply to FEE");
  //  applyFeeButtPtr->MoveResize(   20,  20,  150, 20); 
  
  fWindowPtr->MoveResize(250,250,250,250);
  fWindowPtr->MapSubwindows(); 
  fWindowPtr->MapWindow(); 


}


AliHLTPHOSOnlineDisplayRawDataMenu::~AliHLTPHOSOnlineDisplayRawDataMenu()
{

}

int
AliHLTPHOSOnlineDisplayRawDataMenu::GetStartZ()
{
  return startZInputPtr->GetIntNumber();
}

int
AliHLTPHOSOnlineDisplayRawDataMenu::GetEndZ()
{
  return endZInputPtr->GetIntNumber();
}


int
AliHLTPHOSOnlineDisplayRawDataMenu::GetStartX()
{
  return startXInputPtr->GetIntNumber();
}

int
AliHLTPHOSOnlineDisplayRawDataMenu::GetEndX()
{
   return endXInputPtr->GetIntNumber();
}

int
AliHLTPHOSOnlineDisplayRawDataMenu::GetGain()
{
   return  gainInputPtr->GetIntNumber();
}
