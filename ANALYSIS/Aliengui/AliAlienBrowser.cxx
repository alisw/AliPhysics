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
//                 AliAlienBrowser class
//   The class that deals with the AliEn browser of the GUI
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------


//ROOT
#include "TApplication.h"
#include "TSystem.h"
#include "TObjString.h"

class TGMenu;
class TGToolBar;
class TG3DLine;
class TGStatusBar;
class TGFileDialog;
class TGButton;
class TGIcon;
class TGDockableFrame;
class TGTab;
class TGListTree;

#include "TGFrame.h"
#include "TGListTree.h"

#include "TGrid.h"
#include "TGridResult.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

#include "AliAlienBrowser.h"

ClassImp(AliAlienBrowser)

//___________________________________________________________________________
AliAlienBrowser::AliAlienBrowser(const TGWindow* p, UInt_t w, UInt_t h, TGFrame* frame, const char* classToConnect, EBrowseType type) : TGCanvas(p, w, h), fFrame(frame), fBrowseType(type) {
  // Constructor.
  
  fListTree = new TGListTree(this, kHorizontalFrame);   
  //   fListTree->Associate(this);
  fListTree->Connect("DoubleClicked(TGListTreeItem*,Int_t)", classToConnect,frame, "OnDoubleClick(TGListTreeItem*,Int_t)");
  
  MapSubwindows();
  Resize();
  MapWindow();
}

//___________________________________________________________________________
AliAlienBrowser::~AliAlienBrowser() {
  // Dtcor.
  
  DeleteWindow();
  SetCleanup(kDeepCleanup);
}

//___________________________________________________________________________
void AliAlienBrowser::OnDoubleClick(TGListTreeItem* item, Int_t btn) {
  // OnDoubleClick at the ListTree, depending on Local or Grid Browsing
  
  if(fBrowseType == kLocalBrowse)
    OnDoubleClickLocal(item,btn);
  else if(fBrowseType == kGridBrowse)
    OnDoubleClickGrid(item,btn);
}

//___________________________________________________________________________
void AliAlienBrowser::OnDoubleClickGrid(TGListTreeItem* item, Int_t btn) {
  // OnDoubleClick at the ListTree
  
  if(!gGrid)
    return;
  
  if ((btn!=kButton1) || !item) 
    return;
  
  if((Bool_t)item->GetUserData()) 
    return;
  
  // use UserData to indicate that item was already browsed
  item->SetUserData((void*)1);
  
  TString filename;
  
  TGridResult * result = gGrid->Ls(GetPath(item), "-F");
  
  int i = 0;
  while(result->GetFileName(i)){
    filename = result->GetFileName(i++);
    
    // if the file is a directory
    if(filename.EndsWith("/")){
      if(filename.CompareTo("..") != 0 && filename.CompareTo(".") != 0){
	fListTree->AddItem(item, filename.Remove(filename.Length()-1)); 
      }
    }
  }    
  
  Refresh();
}

//___________________________________________________________________________
void AliAlienBrowser::OnDoubleClickLocal(TGListTreeItem* item, Int_t btn) {
  // Show contents of directory.
  
  if ((btn!=kButton1) || !item || (Bool_t)item->GetUserData()) return;
  
  // use UserData to indicate that item was already browsed
  item->SetUserData((void*)1);
  
  TSystemDirectory dir(item->GetText(),DirName(item));
  
  TList *files = dir.GetListOfFiles();
  
  if (files) {
    TIter next(files);
    TSystemFile *file;
    TString fname;
    
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (file->IsDirectory()) {
	if ((fname!="..") && (fname!=".")) { // skip it
	  fListTree->AddItem(item,fname);
	}
      } 
    }
    fListTree->SortChildren(item);
    delete files;
  }
  Refresh();
}

//___________________________________________________________________________
void AliAlienBrowser::GotoDir(const char* dirToGo) {
  // Goto the given destination dir
  
  fListTree->ClearHighlighted();
  
  TString destDir(dirToGo);
  
  TObjArray * destDirSplit = destDir.Tokenize("/");
  TString dir;
  TGListTreeItem * item = fListTree->GetFirstItem();
  
  if(strcmp(item->GetText(), "/") == 0){
    OnDoubleClick(item, kButton1);
    fListTree->OpenItem(item);
    item = item->GetFirstChild();
  }
  Bool_t found = false;
  
  for(Int_t i=0;i != destDirSplit->GetEntries();i++){
    found = false;
    dir = ((TObjString*)destDirSplit->At(i))->GetString();
    
    do {
      found = dir.CompareTo(item->GetText()) == 0 ? true : false;
      
      if(found){
	
	OnDoubleClick(item, kButton1); // add subdirectories
	
	fListTree->SetSelected(item);
	fListTree->OpenItem(item);
	fListTree->AdjustPosition(item);
	Refresh();
	
      }else{
	item = item->GetNextSibling();
      }
    }while(!found);
    
    item = item->GetFirstChild();  
  }
  
  //   SetHsbPosition(50);
  
  Refresh();
  
  delete destDirSplit;
}

//___________________________________________________________________________
TString AliAlienBrowser::DirName(TGListTreeItem* item) const {
  // Returns an absolute path.
  
  TGListTreeItem* parent;
  TString dirname = item->GetText();
  
  while ((parent=item->GetParent())) {
    dirname = gSystem->ConcatFileName(parent->GetText(),dirname);
    item = parent;
  }
  
  return dirname;
}

//___________________________________________________________________________
void AliAlienBrowser::AddItem(TGListTreeItem* parent, const char* txt) {
  // Add item
  fListTree->AddItem(parent, txt);
  fListTree->SetSelected(fListTree->GetFirstItem());
  Refresh();
}

//___________________________________________________________________________
void AliAlienBrowser::Refresh() const {
  // Refresh the windows
  
  gClient->NeedRedraw(fListTree);
}

//___________________________________________________________________________
const char* AliAlienBrowser::GetPath(TGListTreeItem *item) const {
  // Get the selected item's path. 
  
  //   TGListTreeItem *item = fListTree->GetSelected();
  TString dirName = DirName(item);
  TSystemDirectory dir(item->GetText(),dirName);   
  
  return dirName.Data();
}

//___________________________________________________________________________
const char* AliAlienBrowser::GetPath() {
  // Get the selected item's path. 
  
  return GetPath(fListTree->GetSelected());
}
