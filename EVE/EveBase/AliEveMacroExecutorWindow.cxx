// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMacroExecutorWindow.h"
#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"
#include "AliEveEventManager.h"

#include <TGedEditor.h>
#include <TGLabel.h>
#include <TGListBox.h>
#include <TGTextEntry.h>
#include <Buttons.h>

#include <TPRegexp.h>

class AliEveMEWEntry : public TGTextLBEntry
{
public:
  static void SetFont() { fgDefaultFont = gClient->GetFontPool()->GetFont("-*-lucidatypewriter-*-*-*-*-12-*-*-*-*-*-iso8859-1"); }
};

class AliEveMEWEditor : public TGedEditor
{
protected:
  AliEveMacroExecutorWindow* fMEW;

public:
  AliEveMEWEditor(AliEveMacroExecutorWindow* w) : TGedEditor(0), fMEW(w) {}
  virtual ~AliEveMEWEditor() {}
  virtual void Update(TGedFrame* gframe=0)
  {
    TGedEditor::Update(gframe);
    fMEW->PopulateMacros();
  }
  virtual void Refresh()
  {
    SetModel(fPad, fModel, kButton1Down);
  }
};

//______________________________________________________________________________
// Full description of AliEveMacroExecutorWindow
//

ClassImp(AliEveMacroExecutorWindow)

//______________________________________________________________________________
AliEveMacroExecutorWindow::AliEveMacroExecutorWindow(AliEveMacroExecutor* master) :
  TGMainFrame(gClient->GetRoot()), fM(master),
  fMainFrame(0), fCtrlFrame(0), fListBox(0), fEditor(0),
  fSelectTags(0)
{
  // Constructor.

  fMainFrame = new TGVerticalFrame(this);
  AddFrame(fMainFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

  TGHorizontalFrame *f = 0;
  TGButton* b = 0;
  {
    fCtrlFrame = MkHFrame(fMainFrame);
    fCtrlFrame->Resize(400, 22);

    b = new TGTextButton(fCtrlFrame, "Reload event");
    fCtrlFrame->AddFrame(b);
    b->Connect("Clicked()", "AliEveMacroExecutorWindow", this,
	       "DoReloadEvent()");

    MkLabel(fCtrlFrame, "Select:", 64);
    fSelectTags = new TGTextEntry(f);
    f->AddFrame(fSelectTags, new TGLayoutHints(kLHintsNormal));//|kLHintsExpandX));
    fSelectTags->Connect("ReturnPressed()", "AliEveMacroEditor", this,
			 "DoSelectTags()");
    b = new TGTextButton(fCtrlFrame, "Select");
    fCtrlFrame->AddFrame(b);
    b->Connect("Clicked()", "AliEveMacroExecutorWindow", this,
	       "DoSelectTags()");

    b = new TGTextButton(fCtrlFrame, "Enable all");
    fCtrlFrame->AddFrame(b);
    b->Connect("Clicked()", "AliEveMacroExecutorWindow", this,
	       "DoEnableAll()");

    b = new TGTextButton(fCtrlFrame, "Disable all");
    fCtrlFrame->AddFrame(b);
    b->Connect("Clicked()", "AliEveMacroExecutorWindow", this,
	       "DoDisableAll()");
  }

  fListBox = new TGListBox(fMainFrame);
  fMainFrame->AddFrame(fListBox, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fListBox->Resize(400, 400);
  fListBox->Connect("Selected(Int_t)", "AliEveMacroExecutorWindow", this,
		    "DoMacroSelected(Int_t)");

  fMainFrame->SetEditDisabled(kEditEnable);
  fMainFrame->SetEditable();
  fEditor = new AliEveMEWEditor(this);
  fEditor->SetGlobal(kFALSE);
  fMainFrame->SetEditable(kEditDisable);
  fMainFrame->SetEditable(kFALSE);
  {
    TGFrameElement *el = 0;
    TIter next(fMainFrame->GetList());
    while ((el = (TGFrameElement *) next())) {
      if (el->fFrame == fEditor)
	if (el->fLayout) {
	  el->fLayout->SetLayoutHints(kLHintsExpandX);
	  el->fLayout->SetPadLeft(1); el->fLayout->SetPadRight(1);
	  el->fLayout->SetPadTop(2);  el->fLayout->SetPadBottom(2);
	  break;
	}
    }
  }
  fEditor->Resize(400, 160);
  fEditor->ChangeOptions(fEditor->GetOptions() | kFixedHeight);

  Resize(400, 700);

  SetCleanup(kDeepCleanup);
  Layout();
  MapSubwindows();
  MapWindow();

  gAliEveEvent->Connect("NewEventLoaded()", "AliEveMacroExecutorWindow", this,
			"NewEventLoaded()");
}

AliEveMacroExecutorWindow::~AliEveMacroExecutorWindow()
{
  // Destructor.

  gAliEveEvent->Disconnect("NewEventLoaded()", this);
}

/******************************************************************************/

void AliEveMacroExecutorWindow::PopulateMacros(Bool_t keep_selected)
{
  // Populate list-box (or whatever will replace it) with all macros.
  // prototype: no selection, sorting

  // printf("AliEveMacroExecutorWindow::PopulateMacros()\n");

  AliEveMacro* ex_sel = 0;
  if (keep_selected && fListBox->GetSelected() != -1)
    ex_sel = fBoxContents[fListBox->GetSelected()];

  fListBox->RemoveAll();
  fBoxContents.clear();

  AliEveMEWEntry::SetFont(); 

  TPMERegexp *regexp = 0;
  TString     select = fSelectTags->GetText();
  if ( ! select.IsNull())
  {
    regexp = new TPMERegexp(select, "io");
  }

  TIter next(fM->fMacros);
  AliEveMacro *mac;
  Int_t        id =  0;
  Int_t    sel_id = -1;
  while ((mac = (AliEveMacro*) next()))
  {
    if (regexp && regexp->Match(mac->GetTags()) == 0)
      continue;
    if (mac == ex_sel)
      sel_id = id;

    fListBox->AddEntry(mac->FormForDisplay(), id++);
    fBoxContents.push_back(mac);
  }

  if (sel_id != -1)
    fListBox->Select(sel_id);

  fListBox->MapSubwindows();
  fListBox->Layout();
}

void AliEveMacroExecutorWindow::SetActiveStateOfShownMacros(Bool_t active)
{
  // Set active-state of all shown macros to active.
  // Automatically refreshes the list and editor.

  for (std::vector<AliEveMacro*>::iterator m = fBoxContents.begin(); m != fBoxContents.end(); ++m)
    (*m)->SetActive(active);
  PopulateMacros();
  fEditor->Refresh();
}

/******************************************************************************/

void AliEveMacroExecutorWindow::NewEventLoaded()
{
  // Slot called after a new event has been loaded

  // !!! Once we have exit status from the macro, can update GUI showing this.
}

void AliEveMacroExecutorWindow::DoReloadEvent()
{
  // Slot for reload-event.

  gAliEveEvent->GotoEvent(gAliEveEvent->GetEventId());
}

void AliEveMacroExecutorWindow::DoSelectTags()
{
  // Slot for select-tags.

  PopulateMacros();
}

void AliEveMacroExecutorWindow::DoMacroSelected(Int_t mid)
{
  // Slot for macro-selected.

  fEditor->SetModel(0, fBoxContents[mid], kButton1Down);
}

/******************************************************************************/

TGHorizontalFrame* AliEveMacroExecutorWindow::MkHFrame(TGCompositeFrame* p)
{
  // Make standard horizontal frame.

  if (p == 0)
    p = this;
  TGHorizontalFrame* f = new TGHorizontalFrame(p);
  p->AddFrame(f, new TGLayoutHints(kLHintsNormal|kLHintsExpandX));
  return f;
}

TGLabel* AliEveMacroExecutorWindow::MkLabel(TGCompositeFrame* p, const char* txt, Int_t width,
					    Int_t lo, Int_t ro, Int_t to, Int_t bo)
{
  // Make standard label.

  TGLabel *l = new TGLabel(p, txt);
  p->AddFrame(l, new TGLayoutHints(kLHintsNormal, lo,ro,to,bo));
  l->SetTextJustify(kTextRight);
  l->SetWidth(width);
  l->ChangeOptions(l->GetOptions() | kFixedWidth);
  return l;
}
