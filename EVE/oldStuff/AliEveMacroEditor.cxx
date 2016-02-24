// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMacroEditor.h"
#include "AliEveMacro.h"

#include "TVirtualPad.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGTextEntry.h"

#include "TGComboBox.h"

#include <iostream>

using namespace std;

//______________________________________________________________________________
// GUI editor for AliEveMacro.
//

ClassImp(AliEveMacroEditor)

//______________________________________________________________________________
AliEveMacroEditor::AliEveMacroEditor(const TGWindow *p, Int_t width, Int_t height,
				     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fSources(0),
  fTags(0),
  fMacro(0),
  fFunc(0),
  fArgs(0),
  fActive(0)
{
    cout<<"\n\n\nAliEveMacroEditor constructor called\n\n\n"<<endl;
    
  // Constructor.

  MakeTitle("AliEveMacro");

  TGHorizontalFrame *f = 0;
  //TGLabel           *l = 0;
  Int_t labelW = 48;
  {
    f = MkHFrame();
    MkLabel(f, "Active: ", labelW);

    fActive = new TGCheckButton(f);
    f->AddFrame(fActive); // new TGLayoutHints());
    fActive->Connect("Clicked()", "AliEveMacroEditor", this,
		     "DoActive()");

    MkLabel(f, "Source: ", 56);
    fSources = new TGComboBox(f);
    f->AddFrame(fSources); // new TGLayoutHints());
    fSources->AddEntry("None",      AliEveMacro::kNone);
    fSources->AddEntry("RunLoader", AliEveMacro::kRunLoader);
    fSources->AddEntry("ESD",       AliEveMacro::kESD);
    fSources->AddEntry("ESDfriend", AliEveMacro::kESDfriend);
    fSources->AddEntry("RawReader", AliEveMacro::kRawReader);
    {
      TGListBox* lb = fSources->GetListBox();
      lb->Resize(lb->GetWidth(), 5*16);
    }
    fSources->Resize(92, 20);
    fSources->Connect("Selected(Int_t)", "AliEveMacroEditor", this,
		      "DoSources(Int_t)");

    MkLabel(f, "Tags: ", 40);
    fTags = new TGTextEntry(f);
    f->AddFrame(fTags, new TGLayoutHints(kLHintsNormal|kLHintsExpandX));
    fTags->Connect("TextChanged(const char *)", "AliEveMacroEditor", this,
		    "DoTags()");
  }
  {
    f = MkHFrame();
    MkLabel(f, "Macro: ", labelW);
    fMacro = new TGTextEntry(f);
    f->AddFrame(fMacro, new TGLayoutHints(kLHintsNormal));//|kLHintsExpandX));
    fMacro->Connect("TextChanged(const char *)", "AliEveMacroEditor", this,
		    "DoMacro()");

    MkLabel(f, "Func: ", labelW);
    fFunc = new TGTextEntry(f);
    f->AddFrame(fFunc, new TGLayoutHints(kLHintsNormal|kLHintsExpandX));
    fFunc->Connect("TextChanged(const char *)", "AliEveMacroEditor", this,
		   "DoFunc()");
  }
  {
    f = MkHFrame();
    MkLabel(f, "Args: ", labelW);
    fArgs = new TGTextEntry(f);
    f->AddFrame(fArgs, new TGLayoutHints(kLHintsNormal|kLHintsExpandX));
    fArgs->Connect("TextChanged(const char *)", "AliEveMacroEditor", this,
		   "DoArgs()");
  }
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveMacroEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = static_cast<AliEveMacro*>(obj);

  fSources->Select  (fM->GetSources(), kFALSE);
  fTags   ->SetText (fM->GetTags(),  kFALSE);
  fMacro  ->SetText (fM->GetMacro(), kFALSE);
  fFunc   ->SetText (fM->GetFunc(),  kFALSE);
  fArgs   ->SetText (fM->GetArgs(),  kFALSE);
  fActive ->SetState(fM->fActive ? kButtonDown : kButtonUp);
}

//______________________________________________________________________________
void AliEveMacroEditor::DoSources(Int_t v)
{
   // Slot for Sources.

   fM->SetSources(v);
   Update();
}

//______________________________________________________________________________
void AliEveMacroEditor::DoTags()
{
   // Slot for Tags.

   fM->SetTags(fTags->GetText());
   Update();
}

//______________________________________________________________________________
void AliEveMacroEditor::DoMacro()
{
   // Slot for Macro.

   fM->SetMacro(fMacro->GetText());
   Update();
}

//______________________________________________________________________________
void AliEveMacroEditor::DoFunc()
{
   // Slot for Func.

   fM->SetFunc(fFunc->GetText());
   Update();
}
//______________________________________________________________________________
void AliEveMacroEditor::DoArgs()
{
   // Slot for Args.

   fM->SetArgs(fArgs->GetText());
   Update();
}

//______________________________________________________________________________
void AliEveMacroEditor::DoActive()
{
   // Slot for Active.

   fM->SetActive(fActive->IsOn());
   Update();
}

/******************************************************************************/

TGHorizontalFrame* AliEveMacroEditor::MkHFrame(TGCompositeFrame* p)
{
  // Make standard horizontal frame.

  if (p == 0)
    p = this;
  TGHorizontalFrame* f = new TGHorizontalFrame(p);
  p->AddFrame(f, new TGLayoutHints(kLHintsNormal|kLHintsExpandX));
  return f;
}

TGLabel* AliEveMacroEditor::MkLabel(TGCompositeFrame* p, const char* txt, Int_t width,
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
