// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveGedEditor.h"
#include <TGButton.h>


//==============================================================================
// AliEveGedNameFrame
//==============================================================================

//______________________________________________________________________________
// Full description of AliEveGedNameFrame
//

ClassImp(AliEveGedNameFrame)

//______________________________________________________________________________
AliEveGedNameFrame::AliEveGedNameFrame(const TGWindow *p) :
  TGedFrame(p),
  fB(0)
{
  // Constructor.

  fB = new TGTextButton(this);
  AddFrame(fB, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));
}

//______________________________________________________________________________
void AliEveGedNameFrame::SetModel(TObject* obj)
{
  // Set model object.

  if (obj)
    fB->SetText(Form("%s [%s]", obj->GetName(), obj->ClassName()));
  else
    fB->SetText("No object selected");
}


//==============================================================================
// AliEveGedEditor
//==============================================================================

//______________________________________________________________________________
// Full description of AliEveGedEditor
//

ClassImp(AliEveGedEditor)

//______________________________________________________________________________
AliEveGedEditor::AliEveGedEditor() :
  TEveGedEditor()
{
  // Constructor.

  // Remove old name-frame -- it is created in TGedEditor constructor
  // so virtuals are not active yet.
  fTabContainer->RemoveAll();

  // Replace with a new one.
  TGedFrame* nf = CreateNameFrame(fTabContainer, "Style");
  nf->SetGedEditor(this);
  nf->SetModelClass(0);
  fTabContainer->AddFrame(nf, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
}

//______________________________________________________________________________
TGedFrame* AliEveGedEditor::CreateNameFrame(const TGWindow* parent, const char* /*tab_name*/)
{
  // Create name-frame for a tab.

  return new AliEveGedNameFrame(parent);
}
