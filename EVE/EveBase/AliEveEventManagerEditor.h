// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveEventManagerEditor_H
#define AliEveEventManagerEditor_H

#include "TGedFrame.h"
#include "TGComboBox.h"
#include <TObjString.h>

class AliEveEventManager;
class TEveGValuator;
class TGButton;
class TGCheckButton;
class TGTextButton;
class TGTextView;
class TGNumberEntry;
class TGLabel;

//==============================================================================
// AliEveEventManagerEditor
//==============================================================================

//______________________________________________________________________________
// Short description of AliEveEventManagerEditor
//

class AliEveEventManagerEditor : public TGedFrame
{
public:
   AliEveEventManagerEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
         UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
   virtual ~AliEveEventManagerEditor() {}

   virtual void SetModel(TObject* obj);

   void DoNextEvent();

protected:
   AliEveEventManager  *fM;         // Model object.

   TGTextButton        *fNextEvent; // Load next event
   TGTextView          *fEventInfo; // Text box with event info

private:
   AliEveEventManagerEditor(const AliEveEventManagerEditor&);            // Not implemented
   AliEveEventManagerEditor& operator=(const AliEveEventManagerEditor&); // Not implemented

   ClassDef(AliEveEventManagerEditor, 0); // GUI editor for AliEveEventManager.
};


//==============================================================================
// AliEveEventManagerWindow
//==============================================================================

//______________________________________________________________________________
// Short description of AliEveEventManagerWindow
//

class AliEveEventManagerWindow : public TGMainFrame
{
public:
  AliEveEventManagerWindow(AliEveEventManager* mgr);
  virtual ~AliEveEventManagerWindow();

  void DoFirstEvent();
  void DoPrevEvent();
  void DoNextEvent();
  void DoLastEvent();

  void DoSetEvent();

  void DoRefresh();
  void DoSetAutoLoad();
  void DoSetAutoLoadTime();
  void DoSetTrigSel();

  void Update();

protected:
  AliEveEventManager   *fM;            // Model object.

  TGTextButton         *fFirstEvent;   // Go to first event
  TGTextButton         *fPrevEvent;    // Go to prev event
  TGTextButton         *fNextEvent;    // Go to next event
  TGTextButton         *fLastEvent;    // Go to last event
  TGTextButton         *fRefresh;      // Refresh event-file state

  TGNumberEntry        *fEventId;      // Display/edit current event id
  TGLabel              *fInfoLabel;    // Display last available event id

  TGCheckButton        *fAutoLoad;     // Check-box for automatic loading of events
  TEveGValuator        *fAutoLoadTime; // Time for automatic loading of events

  TGComboBox           *fTrigSel;      // Trigger selection combo box

  TGTextView           *fEventInfo;    // Text box with event info

  TGTextButton* MkTxtButton(TGCompositeFrame* p, const char* txt, Int_t width=0,
			    Int_t lo=0, Int_t ro=0, Int_t to=0, Int_t bo=0);
  TGLabel* MkLabel(TGCompositeFrame* p, const char* txt, Int_t width,
		   Int_t lo=0, Int_t ro=0, Int_t to=2, Int_t bo=0);

private:
  AliEveEventManagerWindow(const AliEveEventManagerWindow&);            // Not implemented
  AliEveEventManagerWindow& operator=(const AliEveEventManagerWindow&); // Not implemented

  ClassDef(AliEveEventManagerWindow, 0); // GUI window for AliEveEventManager.
};

#endif
