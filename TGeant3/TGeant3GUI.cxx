/**************************************************************************
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
*/

/* 
 *  Version: 0
 *  Written by Andreas Morsch
 *  
 * 
 *
 * For questions critics and suggestions to this part of the code
 * contact andreas.morsch@cern.ch
 * 
 **************************************************************************/
#include <stdlib.h>
#include <TROOT.h>
#include <AliRun.h>
#include <AliMC.h>
#include <TGeant3.h>
#include <THIGZ.h>
#include <TApplication.h>
#include <TVirtualX.h>
#include <TGListBox.h>
#include <TGListTree.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGMsgBox.h>
#include <TGMenu.h>
#include <TGCanvas.h>
#include <TGComboBox.h>
#include <TGTab.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>
#include <TGFileDialog.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TEnv.h>
#include <TGPicture.h>

#include "TGeant3GUI.h"

static AliDrawVolume * gCurrentVolume = new AliDrawVolume("NULL");

ClassImp(AliGeant3GeometryGUI)

    AliGeant3GeometryGUI::AliGeant3GeometryGUI()
{
    fPanel  = new AliGuiGeomMain(gClient->GetRoot(), 500, 500);
    fNstack = 0;
    fVolumes = new TClonesArray("AliDrawVolume",100);
    TGeant3 *geant3 = (TGeant3*) gMC;
    fZlq=geant3->Lq();
    fZq=geant3->Q();
    fGclink=geant3->Gclink();
}
void AliGeant3GeometryGUI::Streamer(TBuffer &R__b)
{
;
}

void AliGeant3GeometryGUI::ReadGeometryTree()
{
    char *vname;
    char *namec, *tmp;
    const TGPicture* Folder     = gClient->GetPicture("folder_t.xpm");
    const TGPicture* OpenFolder = gClient->GetPicture("ofolder_t.xpm");
    const TGPicture* Document   = gClient->GetPicture("doc_t.xpm");

    AliDrawVolume  *volume;

    Int_t nst=0;
    Int_t nlevel=1;
    Int_t newlevel=nlevel;

    volume = new AliDrawVolume("ALIC");
    volume->SetIdVolume(((TGeant3*)gMC)->VolId("ALIC"));
    volume->SetIdCopy(0);
    volume->SetItem(NULL);
    (*fVolumes)[0]=volume;

    while(nlevel>nst) {
	for (Int_t i=nst; i<nlevel; i++) 
	{
	    TGListTreeItem *itemi, *item2;
	    Int_t ivol=TMath::Abs(Volume(i)->GetIdVolume());
	    
	    Int_t nch = NChildren(ivol);
	    namec= ((TGeant3*)gMC)->VolName(ivol);
	    if (nch >= 0) {
		printf("\n %s has %d children  \n ", namec,  nch);
	    } else {
		printf("\n %s has %d divisions \n ", namec, -nch);
	    }
	    vname = new char[5];
	    strncpy(vname,namec, 4);
	    vname[4]='\0';
	    Int_t icopy = Volume(i)->GetIdCopy();
	    if (icopy >1) {
		sprintf(namec,"%s*%3dPos",namec,icopy);
	    } else if (icopy <0) {
		sprintf(namec,"%s*%3dDiv",namec,-icopy);
	    }
	    if (i>0) {
		itemi=Volume(i)->GetItem();
	    } else {
		itemi=NULL;
	    }
	    
	    
	    if (nch!=0) {
		item2 = fPanel->AddItem(new AliDrawVolume(vname), 
				    itemi, namec, OpenFolder, Folder);
	    } else {
		item2 = fPanel->AddItem(new AliDrawVolume(vname), 
				    itemi, namec, Document, Document);
	    }
	    
	    if (nch < 0) {
		Int_t icvol=Child(ivol,1);
		namec= ((TGeant3*)gMC)->VolName(-icvol);
		tmp = new char[4];
		strncpy(tmp,(char *) &namec, 4);
		volume = new AliDrawVolume(namec);
		volume->SetIdVolume(-icvol);
		volume->SetIdCopy(nch);
		volume->SetItem(item2);
		(*fVolumes)[newlevel]=volume;
		printf("\n volume %s %d", namec, icvol);
		newlevel++;
	    } else {
		Int_t nnew=0;
		for (Int_t j=0; j<nch; j++) 
		{
		    Int_t icvol=Child(ivol,j+1);
		    icvol = TMath::Abs(icvol);
		    Bool_t inList=kFALSE;
		    for (Int_t k=0; k<nnew; k++) {
			if (icvol==
			    Volume(newlevel-k-1)->GetIdVolume()) 
			{
			    Volume(newlevel-k-1)->AddCopy();
			    inList=kTRUE;
			}
		    }
		    if (!inList) {
			namec=((TGeant3*)gMC)->VolName(icvol);
			tmp = new char[4];
			strncpy(tmp,(char *) &namec, 4);
			volume = new AliDrawVolume(namec);
			volume->SetIdVolume(icvol);
			volume->SetIdCopy(1);
			volume->SetItem(item2);
			(*fVolumes)[newlevel]=volume;
			printf("\n volume %s", namec);
			newlevel++;
			nnew++;
		    }
		}
	    }
	}
	nst=nlevel;
	nlevel=newlevel;
    }
}


Int_t AliGeant3GeometryGUI::NChildren(Int_t idvol)
{
    Int_t jvo = fZlq[fGclink->jvolum-idvol];
    Int_t nin = Int_t(fZq[jvo+3]);
    return nin;
}

Int_t AliGeant3GeometryGUI::Child(Int_t idvol, Int_t idc)
{
    Int_t jvo = fZlq[fGclink->jvolum-idvol];
    Int_t nin=idc;
    Int_t jin = fZlq[jvo-nin];
    Int_t numb =  Int_t (fZq[jin +3]);
    if (numb > 1) {
	return -Int_t(fZq[jin+2]);
    } else {
	return Int_t(fZq[jin+2]);
    }
}


ClassImp(AliDrawVolume)
enum AliDrawParamId {
   P_Theta,
   P_Phi,
   P_Psi,
   P_U,
   P_V,
   P_Uscale,
   P_Vscale,
   P_Shadow,
   P_Hide,
   P_Fill,
   P_Seen,
   P_Clip,
   P_ClipXmin,
   P_ClipXmax,
   P_ClipYmin,
   P_ClipYmax,
   P_ClipZmin,
   P_ClipZmax
};


AliDrawVolume::AliDrawVolume(char* name)
{
    fName   = name;
    fTheta  = 30;
    fPhi    = 30;
    fPsi    = 0;
    fU      = 10;
    fV      = 10;
    fUscale = 0.01;
    fVscale = 0.01;
    fHide=0;
    fShadow=0;
    fFill=1;
    fSeen=1;
    fClip=0;
    fClipXmin=0.;
    fClipXmax=2000.;
    fClipYmin=0.;
    fClipYmax=2000.;
    fClipZmin=0.;
    fClipZmax=2000.;
}

char* AliDrawVolume::Name()
{
    return fName;
}

    
void AliDrawVolume::Streamer(TBuffer &R__b)
{
;
}



void AliDrawVolume::Draw()
{
    gMC->Gsatt(fName,"seen", fSeen);
    
    if (fHide) {
	gMC->Gdopt("hide", "on");
    } else {
	gMC->Gdopt("hide", "off");
    }

    if (fShadow) {
	gMC->Gdopt("shad", "on");
	gMC->Gsatt("*", "fill", fFill);
    } else {
	gMC->Gdopt("shad", "off");
    }

    	gMC->SetClipBox(".");
    if (fClip) {
	gMC->SetClipBox("*", fClipXmin, fClipXmax, fClipYmin, fClipYmax, fClipZmin, fClipZmax);
    } else {
	gMC->SetClipBox(".");
    }
    

    gMC->Gdraw(fName, fTheta, fPhi, fPsi, fU, fV, fUscale, fVscale);
    THIGZ *higz = (THIGZ*)gROOT->GetListOfCanvases()->FindObject("higz");
    if (higz) higz->Update();
}

void AliDrawVolume::DrawSpec()
{
    gMC->Gsatt(fName,"seen", fSeen);
    
    if (fHide) {
	gMC->Gdopt("hide", "on");
    } else {
	gMC->Gdopt("hide", "off");
    }

    if (fShadow) {
	gMC->Gdopt("shad", "on");
	gMC->Gsatt("*", "fill", fFill);
    } else {
	gMC->Gdopt("shad", "off");
    }

    gMC->SetClipBox(".");
    if (fClip) {
	gMC->SetClipBox("*", fClipXmin, fClipXmax, fClipYmin, fClipYmax, fClipZmin, fClipZmax);
    } else {
	gMC->SetClipBox(".");
    }
    

    ((TGeant3*) gMC)->DrawOneSpec(fName);
    THIGZ *higz = (THIGZ*)gROOT->GetListOfCanvases()->FindObject("higz");
    if (higz) higz->Update();
}

void AliDrawVolume::SetParam(Int_t ip, Float_t param)
{
    switch (ip) {
    case P_Theta:
	fTheta=param;
	break;
    case P_Phi:
	fPhi=param;
	break;
    case P_Psi:
	fPsi=param;
	break;
    case P_U:
	fU=param;
	break;
    case P_V:
	fV=param;
	break;
    case P_Uscale:
	fUscale=param;
	break;
    case P_Vscale:
	fVscale=param;
	break;
    case P_Hide:
	fHide=Int_t(param);
	break;
    case P_Shadow:
	fShadow=Int_t(param);
	break;
    case P_Fill:
	fFill=Int_t(param);
	break;
    case P_Seen:
	fSeen=Int_t(param);
	break;
    case P_Clip:
	fClip=Int_t(param);
	break;
    case P_ClipXmin:
	fClipXmin=param;
	break;
    case P_ClipXmax:
	fClipXmax=param;
	break;
    case P_ClipYmin:
	fClipYmin=param;
	break;
    case P_ClipYmax:
	fClipYmax=param;
	break;
    case P_ClipZmin:
	fClipZmin=param;
	break;
    case P_ClipZmax:
	fClipZmax=param;
	break;
    }
}

Float_t  AliDrawVolume::GetParam(Int_t ip)
{
    switch (ip) {
    case P_Theta:
	return fTheta;
	break;
    case P_Phi:
	return fPhi;
	break;
    case P_Psi:
	return fPsi;
	break;
    case P_U:
	return fU;
	break;
    case P_V:
	return fV;
	break;
    case P_Uscale:
	return fUscale;
	break;
    case P_Vscale:
	return fVscale;
	break;
    case P_Hide:
	return Float_t(fHide);
	break;
    case P_Shadow:
	return Float_t(fShadow);
	break;
    case P_Fill:
	return Float_t(fFill);
	break;
    case P_Seen:
	return Float_t(fSeen);
	break;
    case P_Clip:
	return Float_t(fClip);
	break;
    case P_ClipXmin:
	return fClipXmin;
	break;
    case P_ClipXmax:
	return fClipXmax;
	break;
    case P_ClipYmin:
	return fClipYmin;
	break;
    case P_ClipYmax:
	return fClipYmax;
	break;
    case P_ClipZmin:
	return fClipZmin;
	break;
    case P_ClipZmax:
	return fClipZmax;
	break;
    default:
	return 0.;
    }
    return 0.;
}


ClassImp(AliGuiGeomMain)

enum ETestCommandIdentifiers {
   M_FILE_OPEN,
   M_FILE_SAVE,
   M_FILE_SAVEAS,
   M_FILE_EXIT,

   M_TEST_DLG,

   M_HELP_CONTENTS,
   M_HELP_SEARCH,
   M_HELP_ABOUT,


   VId1,
   HId1,
   VId2,
   HId2,

   VSId1,
   HSId1,
   VSId2,
   HSId2
};


Int_t mb_button_id[9] = { kMBYes, kMBNo, kMBOk, kMBApply,
                          kMBRetry, kMBIgnore, kMBCancel,
                          kMBClose, kMBDismiss };

EMsgBoxIcon mb_icon[4] = { kMBIconStop, kMBIconQuestion,
                           kMBIconExclamation, kMBIconAsterisk };

const char *filetypes[] = { "All files",     "*",
                            "ROOT files",    "*.root",
                            "ROOT macros",   "*.C",
                            0,               0 };




TGListTreeItem*  AliGuiGeomMain::
AddItem(TObject * obj, TGListTreeItem *parent, const char* name, const TGPicture *open, const TGPicture *closed)
{
    return fLt->AddItem(parent, name, obj, open, closed);
}

AliGuiGeomMain::AliGuiGeomMain(const TGWindow *p, UInt_t w, UInt_t h)
      : TGMainFrame(p, w, h)
{
    fDialog=0;
    
   // Create test main frame. A TGMainFrame is a top level window.
   // Create menubar and popup menus. The hint objects are used to place
   // and group the different menu widgets with respect to eachother.
   fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                                      0, 0, 1, 1);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

   fMenuFile = new TGPopupMenu(gClient->GetRoot());
   fMenuFile->AddEntry("&Open...", M_FILE_OPEN);
   fMenuFile->AddEntry("&Save", M_FILE_SAVE);
   fMenuFile->AddEntry("S&ave as...", M_FILE_SAVEAS);
   fMenuFile->AddEntry("&Close", -1);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("E&xit", M_FILE_EXIT);

   fMenuFile->DisableEntry(M_FILE_SAVEAS);
   fMenuFile->DisableEntry(M_FILE_OPEN);
   fMenuFile->DisableEntry(M_FILE_SAVE);



   fMenuTest = new TGPopupMenu(this);
   fMenuTest->AddLabel("Draw Control");
   fMenuTest->AddSeparator();
   fMenuTest->AddEntry("&Control Panel ...", M_TEST_DLG);
   fMenuTest->AddSeparator();
   fMenuTest->Associate(this);
   
   
   fMenuHelp = new TGPopupMenu(gClient->GetRoot());
   fMenuHelp->AddEntry("&Contents", M_HELP_CONTENTS);
   fMenuHelp->AddEntry("&Search...", M_HELP_SEARCH);
   fMenuHelp->AddSeparator();
   fMenuHelp->AddEntry("&About", M_HELP_ABOUT);

   fMenuFile->DisableEntry(M_HELP_CONTENTS);
   fMenuFile->DisableEntry(M_HELP_SEARCH);
   fMenuFile->DisableEntry(M_HELP_ABOUT);
   // Menu button messages are handled by the main frame (i.e. "this")
   // ProcessMessage() method.
   fMenuFile->Associate(this);
   fMenuTest->Associate(this);
   fMenuHelp->Associate(this);

   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Draw Control", fMenuTest, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);

   AddFrame(fMenuBar, fMenuBarLayout);
   // Create TreeList
   // Create TGCanvas and a canvas container which uses a tile layout manager

   fCanvasWindow = new TGCanvas(this, 400, 240);
   fLt = new TGListTree(fCanvasWindow->GetViewPort(), 10, 10, kHorizontalFrame,
                        fgWhitePixel);
   fLt->Associate(this);
   fCanvasWindow->SetContainer(fLt);

    
   TGLayoutHints *lo = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   AddFrame(fCanvasWindow, lo);



   SetWindowName("AliRoot Geometry Browser");

   MapSubwindows();

   // we need to use GetDefault...() to initialize the layout algorithm...
   Resize(GetDefaultSize());
   //Resize(400, 200);

   MapWindow();
}

AliGuiGeomMain::~AliGuiGeomMain()
{
   // Delete all created widgets.

   delete fContainer;
   delete fCanvasWindow;

   delete fMenuBarLayout;
   delete fMenuBarItemLayout;
   delete fMenuBarHelpLayout;

   delete fMenuFile;
   delete fMenuTest;
   delete fMenuHelp;
}

void AliGuiGeomMain::Update()
{
    if (fDialog) {
	fDialog->Update();
    }
}

void AliGuiGeomMain::CloseWindow()
{
   // Got close message for this MainFrame. Calls parent CloseWindow()
   // (which destroys the window) and terminate the application.
   // The close message is generated by the window manager when its close
   // window menu item is selected.

   TGMainFrame::CloseWindow();
   gApplication->Terminate(0);
}

Bool_t AliGuiGeomMain::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   // Handle messages send to the AliGuiGeomMain object. E.g. all menu button messages.
    switch (GET_MSG(msg)) {
    case kC_LISTTREE:
	switch (GET_SUBMSG(msg)) {
	case kCT_ITEMCLICK:
	    TGListTreeItem *item;
	    if (parm1 == kButton1) {
		if ((item = fLt->GetSelected())) 
		{
		    gCurrentVolume=((AliDrawVolume *) item->GetUserData());
		    Update();
		}
	    }

	    if (parm1 == kButton2) {
		TGListTreeItem *item;
		if ((item = fLt->GetSelected())) 
		{
		    ((AliDrawVolume *) item->GetUserData())->DrawSpec();
		    gCurrentVolume=((AliDrawVolume *) item->GetUserData());
		    Update();
		}
	    }
	    
	    if (parm1 == kButton3) {
		TGListTreeItem *item;
		if ((item = fLt->GetSelected())) 
		{
		    ((AliDrawVolume *) item->GetUserData())->Draw();
		    gCurrentVolume=((AliDrawVolume *) item->GetUserData());
		    Update();
		}
	    }


	    break;
	case kCT_ITEMDBLCLICK:
	    if (parm1 == kButton1) {
		if (fLt->GetSelected() != 0) {
		    gClient->NeedRedraw(fLt);
		}
	    }
	    break;
	default:
	    break;
	}
	break;
    case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {

            case kCM_BUTTON:
               break;

            case kCM_MENUSELECT:
               break;

            case kCM_MENU:
               switch (parm1) {

                  case M_FILE_OPEN:
                     {
                        TGFileInfo fi;
                        fi.fFileTypes = (char **)filetypes;
                        new TGFileDialog(gClient->GetRoot(), this, kFDOpen,&fi);
                     }
                     break;

                  case M_TEST_DLG:
                     fDialog = new AliGuiGeomDialog
			 (gClient->GetRoot(), this, 400, 200);
                     break;

                  case M_FILE_SAVE:
                     printf("M_FILE_SAVE\n");
                     break;

                  case M_FILE_EXIT:
                     CloseWindow();   // this also terminates theApp
                     break;
                  default:
                     break;
               }
            default:
               break;
         }
      default:
         break;
   }
   return kTRUE;
}

//ClassImp(AliGuiGeomDialog)

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
// Tab2: Seen Option
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
    Float_t param;
//  Update Sliders
    if (fF1) {
	fF1->Update();
    }
//  Seen
    if (fCombo2) {
	param=gCurrentVolume->GetParam(P_Seen);
	fCombo2->Select(Int_t(param)+3);
    }
//  Hide, Shadow, Clip
    if (fChk1) {
	if (Int_t(gCurrentVolume->GetParam(P_Shadow))) {
	    fChk1->SetState(kButtonDown);
	} else {
	    fChk1->SetState(kButtonUp);
	}
    }

    if (fChk2) {
	if (Int_t(gCurrentVolume->GetParam(P_Hide))) {
	    fChk2->SetState(kButtonDown);
	} else {
	    fChk2->SetState(kButtonUp);
	}
    }

    if (fChk3) {
	if (Int_t(gCurrentVolume->GetParam(P_Clip))) {
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
		gCurrentVolume->SetParam(P_ClipXmin,min);
		gCurrentVolume->SetParam(P_ClipXmax,max);
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
		gCurrentVolume->SetParam(P_ClipYmin,min);
		gCurrentVolume->SetParam(P_ClipYmax,max);
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
		gCurrentVolume->SetParam(P_ClipZmin,min);
		gCurrentVolume->SetParam(P_ClipZmax,max);
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
		gCurrentVolume->SetParam(P_Fill, Float_t(parm2));
		gCurrentVolume->Draw();
		break;
	    case 90:
		gCurrentVolume->SetParam(P_Seen, Float_t(parm2-3));
		gCurrentVolume->Draw();
		break;
	    }
	    break;
	case kCM_CHECKBUTTON:
	    switch (parm1) {
	    case 86:
		if (Int_t(gCurrentVolume->GetParam(P_Shadow))) {
		    gCurrentVolume->SetParam(P_Shadow, 0.);
		} else {
		    gCurrentVolume->SetParam(P_Shadow, 1.);
		}
		gCurrentVolume->Draw();
		break;
	    case 87:
		if (Int_t(gCurrentVolume->GetParam(P_Hide))) {
		    gCurrentVolume->SetParam(P_Hide, 0.);
		} else {
		    gCurrentVolume->SetParam(P_Hide, 1.);
		}
		gCurrentVolume->Draw();
		break;
	    case 88:
		if (Int_t(gCurrentVolume->GetParam(P_Clip))) {
		    gCurrentVolume->SetParam(P_Clip, 0.);
		} else {
		    gCurrentVolume->SetParam(P_Clip, 1.);
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

//ClassImp(AliGUISliders)

static Text_t* LabelText[7]  = 
{"Theta  ", "Phi    ", "Psi    ", "U      ", "V      ", "UScale", "VScale"};
static Int_t   IRangeMin[7]  = {    0,     0,     0,    0,    0,   0,   0};
static Int_t   IRangeMax[7]  = {36000, 36000, 36000, 2000, 2000, 10, 10};
static Int_t   DefaultPos[7] = { 3000,  4000,     0, 1000, 1000,   1,   1};

AliGUISliders::AliGUISliders(const TGWindow *p, const TGWindow *main,
                         UInt_t w, UInt_t h) :
    TGCompositeFrame(p, w, h,kVerticalFrame)
{
    ChangeOptions((GetOptions() & ~kHorizontalFrame) | kVerticalFrame);
   //--- layout for buttons: top align, equally expand horizontally
    fBly = new TGLayoutHints(kLHintsTop | kLHintsExpandY, 5, 5, 5, 5);

   //--- layout for the frame: place at bottom, right aligned
    fBfly1 = new TGLayoutHints(kLHintsLeft | kLHintsExpandX );
//
// Frames

   for (Int_t i=0; i<7; i++) {
       Int_t idT=i+1;
       Int_t idS=i+8;       
       fHframe[i] = new TGHorizontalFrame(this, 400, 100, kFixedWidth);
       fTbh[i] = new TGTextBuffer(10);
       fTeh[i] = new TGTextEntry(fHframe[i], fTbh[i],idT);
       char buf[10];
       sprintf(buf, "%6.2f", Float_t(DefaultPos[i])/100);
       fTbh[i]->AddText(0, buf);
       fTeh[i]->Associate(this);
       
       fHslider[i] = new TGHSlider(fHframe[i], 400, kSlider1 | kScaleBoth, idS);
       fHslider[i]->Associate(this);
       fHslider[i]->SetRange(IRangeMin[i], IRangeMax[i]);
       fHslider[i]->SetPosition(DefaultPos[i]);

       fLabel[i] = new TGLabel(fHframe[i], LabelText[i]);
       
       
       fHframe[i]->AddFrame(fLabel[i], fBfly1);
       fHframe[i]->AddFrame(fTeh[i], fBfly1);
       fHframe[i]->AddFrame(fHslider[i], fBfly1);

       AddFrame(fHframe[i], fBly);
   }
}

AliGUISliders::~AliGUISliders()
{
    delete fBfly1; delete fBly;
   // Delete dialog.
    for (Int_t i=1; i<7; i++) {
	delete fHframe[i];
	delete fHslider[i];
	delete fTeh[i];
	delete fTbh[i]; 
    }
}

void AliGUISliders::Update()
{
    char buf[10];
    
    for (Int_t i=0; i<7; i++) {
	Float_t param = gCurrentVolume->GetParam(i);
//
	fHslider[i]->SetPosition(Int_t(param*100.));
	gClient->NeedRedraw(fHslider[i]);
//
	sprintf(buf, "%6.2f", param);
	fTbh[i]->Clear();
	fTbh[i]->AddText(0, buf);
	gClient->NeedRedraw(fTeh[i]);
//
    }

    
}

void AliGUISliders::CloseWindow()
{
   // Called when window is closed via the window manager.

   delete this;
}

Bool_t AliGUISliders::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   // Process slider messages.

   char buf[10];

   switch (GET_MSG(msg)) {
   case kC_TEXTENTRY:
       switch (GET_SUBMSG(msg)) {
       case kTE_TEXTCHANGED:
	   Int_t idT=Int_t(parm1)-1;
	   fHslider[idT]->SetPosition(atof(fTbh[idT]->GetString())*100);
	   gClient->NeedRedraw(fHslider[idT]);
	   gCurrentVolume->SetParam(idT,atof(fTbh[idT]->GetString()));
	   gCurrentVolume->Draw();
       }
       break;
   case kC_HSLIDER:
       switch (GET_SUBMSG(msg)) {
       case kSL_POS:
	   sprintf(buf, "%6.2f", Float_t(parm2)/100);
	   Int_t idS=Int_t(parm1)-8;
	   fTbh[idS]->Clear();
	   fTbh[idS]->AddText(0, buf);
	   gClient->NeedRedraw(fTeh[idS]);
	   gCurrentVolume->SetParam(idS, Float_t(parm2)/100.);
	   gCurrentVolume->Draw();
       }
       break;
   }
   return kTRUE;
}


















