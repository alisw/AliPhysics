/* *************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.8  2000/07/12 08:56:32  fca
Coding convention correction and warning removal

Revision 1.7  2000/06/28 21:27:45  morsch
Most coding rule violations corrected.
Still to do: Split the file (on file per class) ? Avoid the global variables.
Copy constructors and assignment operators (dummy ?)

Revision 1.6  2000/04/14 11:07:46  morsch
Correct volume to medium assignment in case several media are asigned to the
same material.

Revision 1.5  2000/03/20 15:11:03  fca
Mods to make the code compile on HP

Revision 1.4  2000/01/18 16:12:08  morsch
Bug in calculation of number of volume divisions and number of positionings corrected
Browser for Material and Media properties added

Revision 1.3  1999/11/14 14:31:14  fca
Correct small error and remove compilation warnings on HP

Revision 1.2  1999/11/10 16:53:35  fca
The new geometry viewer from A.Morsch

*/

#include "TGButton.h"
#include "TGTab.h"
#include "TGComboBox.h"
#include "TGDoubleSlider.h"

#include "AliGuiGeomDialog.h"
#include "AliGUISliders.h"
#include "AliDrawVolume.h"

AliGuiGeomDialog::AliGuiGeomDialog(const TGWindow *p, const TGWindow *main, UInt_t w,
                       UInt_t h, UInt_t options)
    : TGTransientFrame(p, main, w, h, options)
{
   // Create a dialog window. A dialog window pops up with respect to its
   // "main" window.

   fFrame1 = new TGHorizontalFrame(this, 60, 20, kFixedWidth);

   fOkButton = new TGTextButton(fFrame1, "&Ok", 1);
   fOkButton->Associate(this);
   fCancelButton = new TGTextButton(fFrame1, "&Cancel", 2);
   fCancelButton->Associate(this);

   fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                           2, 2, 2, 2);
   fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1);

   fFrame1->AddFrame(fOkButton, fL1);
   fFrame1->AddFrame(fCancelButton, fL1); 

   fFrame1->Resize(150, fOkButton->GetDefaultHeight());
   AddFrame(fFrame1, fL2);

   //--------- create Tab widget and some composite frames for Tab testing

   fTab = new TGTab(this, 300, 300);
   fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);
//
// Tab1: Sliders
//
   TGCompositeFrame *tf = fTab->AddTab("Draw");
   fF1 = new AliGUISliders(tf, this, 60, 20);
   tf->AddFrame(fF1,fL3);
   
// 
// Tab2: Drawing Options
//
   tf = fTab->AddTab("Options");
   fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                           200, 2, 2, 2);
   fF2 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);

   fF2->AddFrame(fChk1 = new TGCheckButton(fF2, "Shadow", 86), fL3);
   fF2->AddFrame(fChk2 = new TGCheckButton(fF2, "Hide  ", 87), fL3);
   fF2->AddFrame(fChk3 = new TGCheckButton(fF2, "Clip  ", 88), fL3);

   fF2->AddFrame(fLabel1 = new TGLabel(fF2, "Fill"),fL3);
   
   fCombo = new TGComboBox(fF2, 89);
   fF2->AddFrame(fCombo, fL3);

   tf->AddFrame(fF2, fL3);

   int i;
   for (i = 0; i < 8; i++) {
      char tmp[20];

      sprintf(tmp, "%i", i+1);
      fCombo->AddEntry(tmp, i+1);
   }
   fCombo->Select(1);
   fCombo->Resize(50, 20);
   fCombo->Associate(this);

   fChk1->Associate(this);
   fChk2->Associate(this);
   fChk3->Associate(this);
// 
// Tab3: Seen Option
//
   tf = fTab->AddTab("Seen");
   fF3 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
   fF3->AddFrame(fLabel2 = new TGLabel(fF3, "Seen"),fL3);
   fCombo2 = new TGComboBox(fF3, 90);
   fF3->AddFrame(fCombo2, fL3);
   tf->AddFrame(fF3, fL3);

   for (i = 0; i < 4; i++) {
      char tmp[20];

      sprintf(tmp, "%i", i-2);
      fCombo2->AddEntry(tmp, i+1);
   }
   fCombo2->Select(4);
   fCombo2->Resize(50, 20);
   fCombo2->Associate(this);
// 
// Tab4: Clip Box
//
   tf = fTab->AddTab("ClipBox");
   //--- layout for buttons: top align, equally expand horizontally
   fBly = new TGLayoutHints(kLHintsTop | kLHintsExpandY, 5, 5, 5, 5);

   //--- layout for the frame: place at bottom, right aligned
   fBfly1 = new TGLayoutHints(kLHintsLeft | kLHintsExpandX );
//
//  Frames
//
//  Slider1
   fF4 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
       
   fHSframe1 = new TGHorizontalFrame(fF4, 400, 100, kFixedWidth);

   fTbh11 = new TGTextBuffer(10);
   fTeh11 = new TGTextEntry(fHSframe1, fTbh11,1);
   fTbh11->AddText(0, "   0");
   fTeh11->Associate(this);

   fTbh12 = new TGTextBuffer(10);
   fTeh12 = new TGTextEntry(fHSframe1, fTbh12,1);
   fTbh12->AddText(0, "2000");
   fTeh12->Associate(this);
    
   fDslider1 = new TGDoubleHSlider(fHSframe1, 400, kSlider1 | kScaleBoth, 1);
   fDslider1->Associate(this);
   fDslider1->SetRange(-2000, 2000);
   fDslider1->SetPosition(0, 2000);
   
   fSLabel1 = new TGLabel(fHSframe1, "xmin-xmax");

   fHSframe1->AddFrame(fSLabel1, fBfly1);
   fHSframe1->AddFrame(fTeh11, fBfly1);
   fHSframe1->AddFrame(fTeh12, fBfly1);
   fHSframe1->AddFrame(fDslider1, fBfly1);
   
   fF4->AddFrame(fHSframe1, fBly);

//
   fHSframe2 = new TGHorizontalFrame(fF4, 400, 100, kFixedWidth);

   fTbh21 = new TGTextBuffer(10);
   fTeh21 = new TGTextEntry(fHSframe2, fTbh21,1);
   fTbh21->AddText(0, "   0");
   fTeh21->Associate(this);

   fTbh22 = new TGTextBuffer(10);
   fTeh22 = new TGTextEntry(fHSframe2, fTbh22,1);
   fTbh22->AddText(0, "2000");
   fTeh22->Associate(this);

   fDslider2 = new TGDoubleHSlider(fHSframe2, 400, kSlider1 | kScaleBoth, 2);
   fDslider2->Associate(this);
   fDslider2->SetRange(-2000, 2000);
   fDslider2->SetPosition(0, 2000);
   
   fSLabel2 = new TGLabel(fHSframe2, "ymin-ymax");

   fHSframe2->AddFrame(fSLabel2, fBfly1);
   fHSframe2->AddFrame(fTeh21, fBfly1);
   fHSframe2->AddFrame(fTeh22, fBfly1);
   fHSframe2->AddFrame(fDslider2, fBfly1);
   
   fF4->AddFrame(fHSframe2, fBly);

//
   fHSframe3 = new TGHorizontalFrame(fF4, 400, 100, kFixedWidth);

   fTbh31 = new TGTextBuffer(10);
   fTeh31 = new TGTextEntry(fHSframe3, fTbh31,1);
   fTbh31->AddText(0, "   0");
   fTeh31->Associate(this);

   fTbh32 = new TGTextBuffer(10);
   fTeh32 = new TGTextEntry(fHSframe3, fTbh32,1);
   fTbh32->AddText(0, "2000");
   fTeh32->Associate(this);

   fDslider3 = new TGDoubleHSlider(fHSframe3, 400, kSlider1 | kScaleBoth, 3);
   fDslider3->Associate(this);
   fDslider3->SetRange(-2000, 2000);
   fDslider3->SetPosition(0, 2000);
   
   fSLabel3 = new TGLabel(fHSframe3, "zmin-zmax");

   fHSframe3->AddFrame(fSLabel3, fBfly1);
   fHSframe3->AddFrame(fTeh31, fBfly1);
   fHSframe3->AddFrame(fTeh32, fBfly1);
   fHSframe3->AddFrame(fDslider3, fBfly1);
   
   fF4->AddFrame(fHSframe3, fBly);
   tf->AddFrame(fF4, fL3);
//
//
   TGLayoutHints *fL5 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1);
   AddFrame(fTab, fL5);

   MapSubwindows();
   Resize(GetDefaultSize());

   // position relative to the parent's window
   Window_t wdum;
   int ax, ay;
   gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
                          (((TGFrame *) main)->GetWidth() - fWidth) >> 1,
                          (((TGFrame *) main)->GetHeight() - fHeight) >> 1,
                          ax, ay, wdum);
   Move(ax, ay);

   SetWindowName("Dialog");

   MapWindow();
   //gClient->WaitFor(this);    // otherwise canvas contextmenu does not work
}

AliGuiGeomDialog::~AliGuiGeomDialog()
{
   // Delete test dialog widgets.

   delete fOkButton;
   delete fCancelButton;
   delete fFrame1;
   delete fChk1; delete fChk2;
   delete fF1; delete fF2; delete fF3; delete fF4;
   delete fCombo;
   delete fTab;
   delete fL3; delete fL4;
   delete fL1; delete fL2;
   delete fBly; delete fBfly1;
   delete fTbh11; delete fTbh12; delete fTbh21; delete fTbh22; 
   delete fTbh31; delete fTbh32; delete fTeh11; delete fTeh12; 
   delete fTeh21; delete fTeh22; delete fTeh31; delete fTeh32;
   delete fDslider1; delete fDslider2; delete fDslider3;
   delete fSLabel1;  delete fSLabel2;  delete fSLabel3;
}

void AliGuiGeomDialog::Update()
{
// Update widgets
    
    Float_t param;
//  Update Sliders
    if (fF1) {
	fF1->Update();
    }
//  Seen
    if (fCombo2) {
	param=gCurrentVolume->GetParam(kSeen);
	fCombo2->Select(Int_t(param)+3);
    }
//  Hide, Shadow, Clip
    if (fChk1) {
	if (Int_t(gCurrentVolume->GetParam(kShadow))) {
	    fChk1->SetState(kButtonDown);
	} else {
	    fChk1->SetState(kButtonUp);
	}
    }

    if (fChk2) {
	if (Int_t(gCurrentVolume->GetParam(kHide))) {
	    fChk2->SetState(kButtonDown);
	} else {
	    fChk2->SetState(kButtonUp);
	}
    }

    if (fChk3) {
	if (Int_t(gCurrentVolume->GetParam(kClip))) {
	    fChk3->SetState(kButtonDown);
	} else {
	    fChk3->SetState(kButtonUp);
	}
    }
    
}

void AliGuiGeomDialog::CloseWindow()
{
   // Called when window is closed via the window manager.
   delete this;
}

Bool_t AliGuiGeomDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   // Process messages coming from widgets associated with the dialog.
    char buf[10];
    Float_t min,max;
    switch (GET_MSG(msg)) {
    case kC_HSLIDER:
	switch (GET_SUBMSG(msg)) {
	case kSL_POS:
	    switch (Int_t(parm1)) {
	    case 1:
		min=fDslider1->GetMinPosition();
		max=fDslider1->GetMaxPosition();
		sprintf(buf, "%6.2f", min);
		fTbh11->Clear();
		fTbh11->AddText(0, buf);
		sprintf(buf, "%6.2f", max);
		fTbh12->Clear();
		fTbh12->AddText(0, buf);
		gClient->NeedRedraw(fTeh11);
		gClient->NeedRedraw(fTeh12);
		gCurrentVolume->SetParam(kClipXmin,min);
		gCurrentVolume->SetParam(kClipXmax,max);
		break;
	    case 2:
		min=fDslider2->GetMinPosition();
		max=fDslider2->GetMaxPosition();
		sprintf(buf, "%6.2f", min);
		fTbh21->Clear();
		fTbh21->AddText(0, buf);
		sprintf(buf, "%6.2f", max);
		fTbh22->Clear();
		fTbh22->AddText(0, buf);
		gClient->NeedRedraw(fTeh21);
		gClient->NeedRedraw(fTeh22);
		gCurrentVolume->SetParam(kClipYmin,min);
		gCurrentVolume->SetParam(kClipYmax,max);
		break;
	    case 3:
		min=fDslider3->GetMinPosition();
		max=fDslider3->GetMaxPosition();
		sprintf(buf, "%6.2f", min);
		fTbh31->Clear();
		fTbh31->AddText(0, buf);
		sprintf(buf, "%6.2f", max);
		fTbh32->Clear();
		fTbh32->AddText(0, buf);
		gClient->NeedRedraw(fTeh31);
		gClient->NeedRedraw(fTeh32);
		gCurrentVolume->SetParam(kClipZmin,min);
		gCurrentVolume->SetParam(kClipZmax,max);
		break;
	    default:
		break;
	    }
	}
	break;
    case kC_COMMAND:
	switch (GET_SUBMSG(msg)) {
	case kCM_BUTTON:
	    switch(parm1) {
	    case 1:
	    case 2:
		printf("\nTerminating dialog: %s pressed\n",
		       (parm1 == 1) ? "OK" : "Cancel");
		CloseWindow();
		break;
	    }
	    break;
	case kCM_COMBOBOX:
	    switch(parm1) {
	    case 89:
		gCurrentVolume->SetParam(kFill, Float_t(parm2));
		gCurrentVolume->Draw();
		break;
	    case 90:
		gCurrentVolume->SetParam(kSeen, Float_t(parm2-3));
		gCurrentVolume->Draw();
		break;
	    }
	    break;
	case kCM_CHECKBUTTON:
	    switch (parm1) {
	    case 86:
		if (Int_t(gCurrentVolume->GetParam(kShadow))) {
		    gCurrentVolume->SetParam(kShadow, 0.);
		} else {
		    gCurrentVolume->SetParam(kShadow, 1.);
		}
		gCurrentVolume->Draw();
		break;
	    case 87:
		if (Int_t(gCurrentVolume->GetParam(kHide))) {
		    gCurrentVolume->SetParam(kHide, 0.);
		} else {
		    gCurrentVolume->SetParam(kHide, 1.);
		}
		gCurrentVolume->Draw();
		break;
	    case 88:
		if (Int_t(gCurrentVolume->GetParam(kClip))) {
		    gCurrentVolume->SetParam(kClip, 0.);
		} else {
		    gCurrentVolume->SetParam(kClip, 1.);
		}
		gCurrentVolume->Draw();
		break;

	    default:
		break;
	    }
	    break;
	case kCM_TAB:
	    break;
	default:
	    break;
	}
	break;
    default:
	break;
    }
    return kTRUE;
}

