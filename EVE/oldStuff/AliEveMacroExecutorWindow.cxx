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

class AliEveMEWListBox : public TGListBox
{
public:
  AliEveMEWListBox(const TGWindow* p = 0, Int_t id = -1) : TGListBox(p, id)
  {
    if (gfGC == 0)
    {
      const TGGC& old = TGTextLBEntry::GetDefaultGC();

      gfFont = gClient->GetFontPool()->GetFont("-*-lucidatypewriter-*-*-*-*-12-*-*-*-*-*-iso8859-1");
      gfGC   = gClient->GetGCPool()->GetGC((GCValues_t*) old.GetAttributes(), kTRUE);
      gfGC->SetFont(gVirtualX->GetFontHandle(gfFont->GetFontStruct()));
    }
  }

  using TGListBox::AddEntry;
  virtual void AddEntry(const char* s, Int_t id)
  {
    static const Pixel_t gkBackground[] = { 0x00ffffff, 0xf5f7f8 };

    TGTextLBEntry *lbe = new TGTextLBEntry(fLbc, new TGString(s), id,
					   gfGC->GetGC(), gfFont->GetFontStruct());
    fItemVsize = TMath::Max(fItemVsize, lbe->GetDefaultHeight());
    fLbc->AddEntry(lbe, new TGLayoutHints(kLHintsExpandX | kLHintsTop));
    // Need to set it here as the above line sets it to white (for some strange reason).
    lbe->SetBackgroundColor(gkBackground[id%2]);
  }

protected:
  static TGGC    *gfGC;
  static TGFont  *gfFont;

private:
  AliEveMEWListBox(const AliEveMEWListBox&);            // Not implemented
  AliEveMEWListBox& operator=(const AliEveMEWListBox&); // Not implemented
};

TGGC   *AliEveMEWListBox::gfGC   = 0;
TGFont *AliEveMEWListBox::gfFont = 0;


class AliEveMEWEditor : public TGedEditor
{
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
protected:
  AliEveMacroExecutorWindow* fMEW;
private:
  AliEveMEWEditor(const AliEveMEWEditor&);            // Not implemented
  AliEveMEWEditor& operator=(const AliEveMEWEditor&); // Not implemented
};

//______________________________________________________________________________
// Full description of AliEveMacroExecutorWindow
//

ClassImp(AliEveMacroExecutorWindow)

//______________________________________________________________________________
AliEveMacroExecutorWindow::AliEveMacroExecutorWindow(AliEveMacroExecutor* master) :
  TGMainFrame(gClient->GetRoot()), fM(master),
  fMainFrame(0), fCtrlFrame(0), fListBox(0), fEditor(0),
  fSelectTags(0),
  fBoxContents()
{
  // Constructor.

  fMainFrame = new TGVerticalFrame(this);
  AddFrame(fMainFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

  // TGHorizontalFrame *f = 0;
  TGButton          *b = 0;
  {
    fCtrlFrame = MkHFrame(fMainFrame);
    fCtrlFrame->Resize(400, 22);

    b = new TGTextButton(fCtrlFrame, "Reload event");
    fCtrlFrame->AddFrame(b);
    b->Connect("Clicked()", "AliEveMacroExecutorWindow", this,
	       "DoReloadEvent()");

    MkLabel(fCtrlFrame, "Select: ", 48);
    fSelectTags = new TGTextEntry(fCtrlFrame);
    fCtrlFrame->AddFrame(fSelectTags, new TGLayoutHints(kLHintsNormal));//|kLHintsExpandX));
    fSelectTags->Connect("ReturnPressed()", "AliEveMacroExecutorWindow", this,
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

  fListBox = new AliEveMEWListBox(fMainFrame);
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
  fEditor->Resize(400, 150);
  fEditor->ChangeOptions(fEditor->GetOptions() | kFixedHeight);

  Resize(400, 700);

  //SetCleanup(kDeepCleanup);
  Layout();
  MapSubwindows();
  MapWindow();

  AliEveEventManager::Instance()->Connect("NewEventLoaded()", "AliEveMacroExecutorWindow", this,
			"NewEventLoaded()");
}

AliEveMacroExecutorWindow::~AliEveMacroExecutorWindow()
{
  // Destructor.
	fBoxContents.clear();
//    AliEveEventManager::Instance()->Disconnect("NewEventLoaded");/*()", this);*/
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

  Int_t sbar_pos = fListBox->GetVScrollbar()->GetPosition();

  fListBox->RemoveAll();
  fBoxContents.clear();

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

  fListBox->Layout();
  fListBox->GetVScrollbar()->SetPosition(sbar_pos);
  fListBox->MapSubwindows();
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

  AliEveEventManager::Instance()->GotoEvent(AliEveEventManager::Instance()->GetEventId());
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
