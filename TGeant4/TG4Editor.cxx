// $Id$
// Category: interfaces
//
// Author: D. Adamova
//===================================================
//
//--------TG4Editor.cxx--------------------------//
//------- A service class for GUI--------------//
//
//=================================================

#include "TG4Editor.h"

#include <TGButton.h>
#include <TGTextEdit.h>
#include <TGText.h>


ClassImp(TG4Editor)

TG4Editor::TG4Editor(const TGWindow* main, UInt_t w, UInt_t h) :
    TGTransientFrame(gClient->GetRoot(), main, w, h)
{
   // Create an editor 

   fEdit = new TGTextEdit(this, w, h, kSunkenFrame | kDoubleBorder);
   fL1 = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 3, 3, 3, 3);
   AddFrame(fEdit, fL1);

   fOK = new TGTextButton(this, "  &OK  ");
   fL2 = new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 5, 5);
   AddFrame(fOK, fL2);

   SetTitle();

   MapSubwindows();

   Resize(GetDefaultSize());

   // position relative to the parent's window
   Window_t wdum;
   int ax, ay;
   gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
                          ((TGFrame *) main)->GetWidth() - (fWidth >> 1),
                          (((TGFrame *) main)->GetHeight() - fHeight) >> 1,
                          ax, ay, wdum);
   Move(ax, ay);
   SetWMPosition(ax, ay);
}

TG4Editor::~TG4Editor()
{
   // Delete editor accessories

   delete fEdit;
   delete fOK;
   delete fL1;
   delete fL2;
}

void TG4Editor::SetTitle()
{
   // Set title in editor window.

   TGText* txt = GetEditor()->GetText();
   Bool_t untitled = !strlen(txt->GetFileName()) ? kTRUE : kFALSE;

   char title[256];
   //if (untitled)
      sprintf(title, "Status Report");
   //else
      //sprintf(title, "Editor - %s", txt->GetFileName());

   SetWindowName(title);
   SetIconName(title);
}

void TG4Editor::Popup()
{
   // Show editor.

   MapWindow();
}

void TG4Editor::LoadBuffer(const char* buffer)
{
   // Load a text buffer in the editor.

   fEdit->LoadBuffer(buffer);
}

void TG4Editor::LoadFile(const char* file)
{
   // Load a file in the editor.

   fEdit->LoadFile(file);
}

void TG4Editor::CloseWindow()
{
   // Called when closed via window manager action.

   delete this;
}

Bool_t TG4Editor::ProcessMessage(Long_t msg, Long_t, Long_t)
{
   // Process Help Menu.
   

   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               // Only one button and one action...
               delete this;
               break;
            default:
               break;
         }
         break;
      case kC_TEXTVIEW:
         switch (GET_SUBMSG(msg)) {
            case kTXT_CLOSE:
               // close window
               delete this;
               break;
            case kTXT_OPEN:
               SetTitle();
               break;
            case kTXT_SAVE:
               SetTitle();
               break;
            default:
               break;
         }
      default:
         break;
   }

   return kTRUE;
}

