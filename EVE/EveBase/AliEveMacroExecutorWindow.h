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
class TGListBox;
class TGedEditor;

//______________________________________________________________________________
// Short description of AliEveMacroExecutorWindow
//

class AliEveMacroExecutorWindow : public TGMainFrame
{
public:
  AliEveMacroExecutorWindow(AliEveMacroExecutor* master);
  virtual ~AliEveMacroExecutorWindow();

  void PopulateMacros(Bool_t keep_selected=kTRUE);

  void NewEventLoaded();

  void DoReloadEvent();
  void DoMacroSelected(Int_t mid);

protected:
  AliEveMacroExecutor *fM;

  TGCompositeFrame *fMainFrame;
  TGCompositeFrame *fCtrlFrame;
  TGListBox        *fListBox;
  TGedEditor       *fEditor;

  std::vector<AliEveMacro*> fBoxContents;

private:
  AliEveMacroExecutorWindow(const AliEveMacroExecutorWindow&);            // Not implemented
  AliEveMacroExecutorWindow& operator=(const AliEveMacroExecutorWindow&); // Not implemented

  ClassDef(AliEveMacroExecutorWindow, 0); // Short description.
};

#endif
