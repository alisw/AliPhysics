 // $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveHLTEventManagerEditor_H
#define AliEveHLTEventManagerEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGTextButton;
class TGNumberEntry;
class TGColorSelect;
class TGPictureButton;
class TGComboBox;
class AliEveHLTEventManager;

class AliEveHLTEventManagerEditor : public TGedFrame
{
public:
  AliEveHLTEventManagerEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, 
			   UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveHLTEventManagerEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void ConnectToHLT();
  void NextEvent();
  void EventLoop();

  void NavigateBack();
  void NavigateFwd();
  void SetTriggerString(int id);
  void WriteBlockListToFile();
  void PrintScreens();
  void PollEvents();

protected:

  AliEveHLTEventManager  *fM; // Model object.
  
  TGTextButton     *fButtonConnect; // Button to connect to HOMER.
  TGTextButton     *fButtonWriteToFile; // Button to write block list to file
  TGTextButton     *fButtonNextEvent; // Button to call next Even
  TGTextButton     *fButtonNavigateBack; // Button to navigate back
  TGTextButton     *fButtonNavigateFwd;  // Button to navigate fwd
  TGTextButton     *fButtonPrintScreens;  // Button to print viewers
  //TGComboBox       *fBoxTriggerSelector; // Drop down menu to select trigger criteria.
  TGTextButton     *fButtonEventLoopText; //Text button to start / stop event loop.
  TGTextButton    *fButtonUpdateEvents;
  //TGComboBox       *fBoxEventLoopSpeed; // Drop down menu to set the speed of the loop.
  TGPictureButton  *fButtonEventLoop; // Picture button to start/stop event loop, HLT LOGO.
  

private:
  AliEveHLTEventManagerEditor(const AliEveHLTEventManagerEditor&);            // Not implemented
  AliEveHLTEventManagerEditor& operator=(const AliEveHLTEventManagerEditor&); // Not implemented

  Bool_t fEventLoopStarted;
  Bool_t fBufferLoopStarted;

  ClassDef(AliEveHLTEventManagerEditor, 0); // Editor for AliEveHLTEventManager
};

#endif
