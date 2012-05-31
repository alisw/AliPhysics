// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMacroExecutorWindow_H
#define AliEveMacroExecutorWindow_H

#include "TGFrame.h"

#include <vector>

class AliEveMacroExecutor;
class AliEveMacro;
class AliEveMEWEditor;

class TGLabel;
class TGListBox;
class TGTextEntry;

//______________________________________________________________________________
// Short description of AliEveMacroExecutorWindow
//

class AliEveMacroExecutorWindow : public TGMainFrame
{
public:
  AliEveMacroExecutorWindow(AliEveMacroExecutor* master);
  virtual ~AliEveMacroExecutorWindow();

  void PopulateMacros(Bool_t keep_selected=kTRUE);

  void SetActiveStateOfShownMacros(Bool_t active);

  void NewEventLoaded();

  void DoEnableAll()  { SetActiveStateOfShownMacros(kTRUE);  }
  void DoDisableAll() { SetActiveStateOfShownMacros(kFALSE); }
  void DoReloadEvent();
  void DoSelectTags();
  void DoMacroSelected(Int_t mid);

protected:
  AliEveMacroExecutor *fM;

  TGCompositeFrame *fMainFrame;
  TGCompositeFrame *fCtrlFrame;
  TGListBox        *fListBox;
  AliEveMEWEditor  *fEditor;

  TGTextEntry      *fSelectTags;

  std::vector<AliEveMacro*> fBoxContents;

  TGHorizontalFrame* MkHFrame(TGCompositeFrame* p=0);
  TGLabel*           MkLabel (TGCompositeFrame* p, const char* txt, Int_t width,
			      Int_t lo=0, Int_t ro=0, Int_t to=2, Int_t bo=0);

private:
  AliEveMacroExecutorWindow(const AliEveMacroExecutorWindow&);            // Not implemented
  AliEveMacroExecutorWindow& operator=(const AliEveMacroExecutorWindow&); // Not implemented

  ClassDef(AliEveMacroExecutorWindow, 0); // Short description.
};

#endif
