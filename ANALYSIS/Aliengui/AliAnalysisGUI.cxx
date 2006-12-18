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

/* $Id$ */

//-----------------------------------------------------------------
//           AliAnalysisGUI class
//   The class that deals with the analysis GUI
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

//ROOT
#include "TApplication.h"
#include "TSystem.h"

#include "TGTab.h"
#include "TGDockableFrame.h"
#include "TGFrame.h"
#include "TGMenu.h"
#include "TG3DLine.h"
#include "TGToolBar.h"
#include "TGStatusBar.h"
#include "TGIcon.h"
#include "TGFileDialog.h"

class TGButton;

#include "TGrid.h"

#define GUIDEBUG 1

#ifdef GUIDEBUG
using namespace std;
#endif


enum EAliEnViewerCommands {
   kAliEnConnect,
   kAliEnHome,
   kOpen,
   kSave,
   kExit
};

//GUI
#include "AliAlienBrowser.h"
#include "AliFileListFrame.h"
#include "AliLoginFrame.h"
#include "AliPackageFrame.h"

#include "AliAnalysisGUI.h"

ClassImp(AliAnalysisGUI)

//___________________________________________________________________________
AliAnalysisGUI::AliAnalysisGUI(const TGWindow *p, UInt_t w, UInt_t h) : 
  TGMainFrame(p,w,h), 
  fHFrame1(0), fVFrame1(0), fVFrame2(0),
  fMenuDock(0), fMenuFile(0), fMenuBar(0),
  fToolBar(0), fTab(0),
  fMenuBarLayout(0), fMenuBarItemLayout(0),
  fH3DLine(0), fCanvas2(0), fStatusBar(0),
  fAliEnBrowser(0), fFileListFrame(0),
  fLogInFrame(0), fTagFrame(0),
  fTagAnalysisFrame(0), fPackageFrame(0),
  fSelectorFrame(0), fIcon(0),
  fRightIconPicture(0), fRightIcon(0),
  fIsConnected(kFALSE), fAlien(0) {
  // AliAnalysisGUI Constructor
  
  SetWindowName("AliEn");
  
  // Create all the Frames, MenuBar and ToolBar
  fVFrame1 = new TGVerticalFrame(this, 600, 600);
  
  AddMenuBar();
  AddToolBar();
  
  AddFrame(fVFrame1, new TGLayoutHints(kLHintsTop));
  
  fTab = new TGTab(this, 900, 300);
  //   fTab->Connect("Selected(Int_t)", "TestDialog", this, "DoTab(Int_t)");
  
  AddFrame(fTab);
  
  //_____________________________________//
  //_________File Catalogue TAB__________//
  //_____________________________________//
  TGCompositeFrame *tf = fTab->AddTab("File Catalogue");   
  TGCompositeFrame *fF1 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  
  fHFrame1 = new TGHorizontalFrame(fF1, 500, 500);
  fHFrame1->SetCleanup(kDeepCleanup);
  
  fAliEnBrowser = new AliAlienBrowser(fHFrame1, 200, 250, this, "AliAnalysisGUI", kGridBrowse);
  fHFrame1->AddFrame(fAliEnBrowser, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));         
  
  fVFrame2 = new TGVerticalFrame(fHFrame1, 500, 500);
  fVFrame2->SetCleanup(kDeepCleanup);
  fHFrame1->AddFrame(fVFrame2, new TGLayoutHints(kLHintsRight | kLHintsExpandY));                 
  
  fFileListFrame = new AliFileListFrame(fVFrame2, 250, 200);   
  
  fVFrame2->AddFrame(fFileListFrame,new TGLayoutHints(kLHintsExpandY));
  
  fF1->AddFrame(fHFrame1,new TGLayoutHints(kLHintsNormal));      
  tf->AddFrame(fF1, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
  
  
  //_____________________________________//
  //_________Packages TAB__________//
  //_____________________________________//
  tf = fTab->AddTab("Packages");
  
  fPackageFrame = new AliPackageFrame(tf, 250, 300, this);
  tf->AddFrame(fPackageFrame, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));
  
  //_____________________________________//
  //______________TAG TAB________________//
  //_____________________________________//
  tf = fTab->AddTab("Event Tags");
  
  //TGCompositeFrame *fF2 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
  
  fTagAnalysisFrame = new AliTagAnalysisFrame(tf, 250, 200, this);
  
  tf->AddFrame(fTagAnalysisFrame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));
  
  //_____________________________________//
  //___________SELECTOR TAB______________//
  //_____________________________________//
  tf = fTab->AddTab("Analysis");
  
  fSelectorFrame = new AliSelectorFrame(tf, 250, 300, this, fTagAnalysisFrame);
  
  tf->AddFrame(fSelectorFrame, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));
  
  
  // Status Bar
  AddStatusBar();

  MapSubwindows();
  Resize();
  MapWindow();
}

//___________________________________________________________________________
AliAnalysisGUI::~AliAnalysisGUI() {
  // AliAnalysisGUI Destructor
  
  Cleanup();
  
  delete fMenuDock;
  delete fMenuFile;
  delete fMenuBar;
  delete fToolBar;
  delete fH3DLine;
  delete fHFrame1;
  delete fCanvas2;
  delete fFileListFrame;
  delete fIcon;
  delete fMenuBarLayout;
  delete fMenuBarItemLayout;
}

//___________________________________________________________________________
void AliAnalysisGUI::CloseWindow() {
  // Got close message for this MainFrame. Terminates the application.
  
  gApplication->Terminate();
}

//___________________________________________________________________________
void AliAnalysisGUI::AddMenuBar() {
  // Create the MenuBar
  
  fMenuDock = new TGDockableFrame(this);
  AddFrame(fMenuDock, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));
  fMenuDock->SetWindowName("Menu");
  
  fMenuFile = new TGPopupMenu(gClient->GetRoot());
  fMenuFile->AddEntry("&Log In...", kMFILELOGIN);
  fMenuFile->AddEntry("&Open...", kMFILEOPEN);
  fMenuFile->AddEntry("S&ave as...", kMFILESAVEAS);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("Tag...", kMFILETAG);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("E&xit", kMFILEEXIT);
  
  fMenuBar = new TGMenuBar(fMenuDock, 1, 1, kHorizontalFrame);
  fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
  
  fMenuDock->AddFrame(fMenuBar, fMenuBarLayout);
  
  //   AddFrame(fMenuDock, fMenuBarLayout);
  fVFrame1->AddFrame(fMenuDock, fMenuBarLayout);
  
  
  
  //   fMenuDock->Connect("Undocked()", "AliAnalysisGUI", this, "HandleMenu(=M_VIEW_UNDOCK)");
  
  fMenuFile->Connect("Activated(Int_t)", "AliAnalysisGUI", this, "HandleMenu(Int_t)");
  
  MapSubwindows();
  Resize();
  MapWindow();
}

//___________________________________________________________________________
void AliAnalysisGUI::AddToolBar() {
  // Create the ToolBar
  
  // toolbar icon files
  const char *xpmtoolbar[] = {
    "connect.xpm",
    "",
    "home_t.xpm",
    "",
    "fileopen.xpm",
    "filesaveas.xpm",
    "",
    "ed_quit.png",
    "",
  };
  
  ToolBarData_t tbdata[] = {
    { "", "Quick Connect",        kFALSE, kAliEnConnect,    0 },
    { "", 0,                      kFALSE, -1,               0 },
    { "", "Home Directory",       kFALSE, kAliEnHome,       0 },
    { "", 0,                      kFALSE, -1,               0 },
    { "", "Open",                 kFALSE, kOpen,            0 },
    { "", "Save",                 kFALSE, kSave,            0 },
    { "", 0,                      kFALSE, -1,               0 },
    { "", "Exit",                 kFALSE, kExit,            0 },
    { 0,  0,                      0,      0,                0 }
  };
  
  // toolbar button separator
  int separator = 5;
  
  // number of icons
  const int knumIcons= 8;
  
  // creation of a toolbar object as a child of main frame
  fToolBar = new TGToolBar(this, 640, 80);
  for (int i = 0; i < knumIcons; i++) {
    // filling the ToolBarData_t with information
    tbdata[i].fPixmap = xpmtoolbar[i];
    
    if (strlen(xpmtoolbar[i]) == 0) {
      separator = 5;
      continue;
    }
    fToolBar->AddButton(this, &tbdata[i], separator);
    separator = 0;
  }
  
  // adding the tool bar to the main frame
  //   AddFrame(fToolBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
  fVFrame1->AddFrame(fToolBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
  // adding a horizontal line as a separator
  fH3DLine = new TGHorizontal3DLine(this);
  //   AddFrame(fH3DLine, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
  
  fVFrame1->AddFrame(fH3DLine, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
  
  fToolBar->Connect("Clicked(Int_t)", "AliAnalysisGUI", this, "HandleToolBar(Int_t)");
  
  MapSubwindows();
  Resize();
  MapWindow();
}

//___________________________________________________________________________
void AliAnalysisGUI::AddStatusBar() {
  // Status Bar
  Int_t parts[] = {40, 30, 30};
  fStatusBar = new TGStatusBar(this, 50, 10, kHorizontalFrame);
  fStatusBar->SetParts(parts, 3);
  
  AddFrame(fStatusBar, new TGLayoutHints(kLHintsBottom|kLHintsExpandX, 0,0,0,0));
   
  UserGroup_t * fUserGroup = gSystem->GetUserInfo();
  char line [100];
  
  sprintf(line,"User : %s - %s", fUserGroup->fRealName.Data(), fUserGroup->fGroup.Data());
  fStatusBar->SetText(line, 1);
  
  TGCompositeFrame *leftpart = fStatusBar->GetBarPart(2);
  
  if(!fIsConnected){
    fStatusBar->SetText("      Disconnected", 2);
    fRightIconPicture = (TGPicture *)fClient->GetPicture("proof_disconnected.xpm");
  }
  else{
    fStatusBar->SetText("      Connected", 2);
    fRightIconPicture = (TGPicture *)fClient->GetPicture("monitor01.xpm");
  }
  
  fRightIcon = new TGIcon(leftpart, fRightIconPicture,fRightIconPicture->GetWidth(),fRightIconPicture->GetHeight());
  leftpart->AddFrame(fRightIcon, new TGLayoutHints(kLHintsLeft, 2, 0, 0, 0)); 
}

//___________________________________________________________________________
Bool_t AliAnalysisGUI::LogIn(const char * server, const char *username) {
  // Log in to AliEn
  
  //   fAlien = TGrid::Connect(server, username); 
  fAlien = TGrid::Connect(server); 
  
  fIsConnected = gGrid;
  if(fIsConnected){
    fAliEnBrowser->AddItem(0, "/");
    fStatusBar->SetText("      Connected", 2);
    ChangeRightLogo("proof_connected.xpm");
  }
  
  return fIsConnected;
}

//___________________________________________________________________________
void AliAnalysisGUI::ChangeRightLogo(const char *name) {
  // Change the right logo (used for animation).
  
  fClient->FreePicture(fRightIconPicture);
  fRightIconPicture = (TGPicture *)fClient->GetPicture(name);
  fRightIcon->SetPicture(fRightIconPicture);
}

//___________________________________________________________________________
void AliAnalysisGUI::OnDoubleClick(TGListTreeItem* item, Int_t btn) {
  // OnDoubleClick at the ListTree
  
  fAliEnBrowser->OnDoubleClick(item, btn);   
  fFileListFrame->SetQueryPath(fAliEnBrowser->GetPath());
}

//___________________________________________________________________________
void AliAnalysisGUI::HandleMenu(Int_t id) {
  // Handle menu items.
  
  const char *gAlifiletypes[] = { 
    "ROOT files",    "*.root",
    "ROOT macros",   "*.C",
    "Text files",    "*.[tT][xX][tT]",
    "All files",     "*",
    0,               0 
  };

  switch (id) {
  case kMFILELOGIN:
    fLogInFrame = new AliLoginFrame(gClient->GetRoot(), this, 400, 200);
    break;
  case kMFILEOPEN: {
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = gAlifiletypes;
    fi.fIniDir    = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    printf("Open file: %s (dir: %s)\n", fi.fFilename, fi.fIniDir);
    dir = fi.fIniDir;	 
  }
    break;
    
  case kMFILESAVEAS:
    
    break;
  case kMFILETAG:
    //	new TagFrame(gClient->GetRoot(), this, 400, 200, kHorizontalFrame, -1);
    
    break;	
  case kMFILEEXIT:
    CloseWindow();
    break;
  }
}

//___________________________________________________________________________
void AliAnalysisGUI::HandleToolBar(Int_t id) {
  // Handle menu items.
  
  switch (id) {
  case kAliEnConnect:
    
    LogIn("alien://", "");
    break;
  case kAliEnHome:
    if(gGrid){
      fAliEnBrowser->GotoDir(gGrid->GetHomeDirectory());       
      fFileListFrame->SetQueryPath(fAliEnBrowser->GetPath()); 
    }
    break;
  case kExit:
    CloseWindow();
    break;
  }
}
