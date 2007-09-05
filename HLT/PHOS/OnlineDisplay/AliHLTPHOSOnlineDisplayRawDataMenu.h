#ifndef ALIHLTPHOSONLINEDISPLAYRAWDATAMENU_H
#define ALIHLTPHOSONLINEDISPLAYRAWDATAMENU_H

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

//class TGWindow;
//class TGFrame;
class TGMainFrame;
class TGGroupFrame;
class TGNumberEntry;
class TGLabel;
class AliHLTPHOSOnlineDisplayNumberEntry;
class AliHLTPHOSOnlineDisplayButton;
class AliHLTPHOSOnlineDispplayEventTab;
class AliHLTPHOSOnlineDisplay;

class  AliHLTPHOSOnlineDisplayRawDataMenu
{
public:
  int GetStartZ();
  int GetEndZ();
  int GetStartX();
  int GetEndX();
  int GetGain();

  //  AliHLTPHOSOnlineDisplayRawDataMenu();
  AliHLTPHOSOnlineDisplayRawDataMenu(AliHLTPHOSOnlineDisplay *onlineDisplayPtr);

  virtual ~AliHLTPHOSOnlineDisplayRawDataMenu();
  
  //  TGWindow *fWindowPtr;
  TGMainFrame *fWindowPtr;
  //  TGGroupFrame  *fFramePtr;
  AliHLTPHOSOnlineDisplayNumberEntry   *startZInputPtr;
  AliHLTPHOSOnlineDisplayNumberEntry   *endZInputPtr;
  AliHLTPHOSOnlineDisplayNumberEntry   *startXInputPtr;
  AliHLTPHOSOnlineDisplayNumberEntry   *endXInputPtr;

  AliHLTPHOSOnlineDisplayNumberEntry   *gainInputPtr; 

  TGLabel         *zLabelPtr;
  TGLabel         *xLabelPtr;
  TGLabel         *startLabelPtr;
  TGLabel         *endLabelPtr;

  AliHLTPHOSOnlineDisplayButton *fGetDataButtonPtr;
  TGLabel         *gainLabelPtr;
  //  AliHLTPHOSOnlineDispplayEventTab *fEventPtr;

private:
  AliHLTPHOSOnlineDisplayRawDataMenu();
  AliHLTPHOSOnlineDisplay *fOnlineDisplayPtr;

};

#endif
