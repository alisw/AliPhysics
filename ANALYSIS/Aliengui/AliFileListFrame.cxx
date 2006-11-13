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
//           AliFileListFrame class
//   The class that deals with the file list frame of the GUI
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "TObjString.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TGMsgBox.h"
#include "TGTableLayout.h"
#include "TGListView.h"

#include "TFile.h"

#include "TGrid.h"
#include "TGridResult.h"

#include "AliFileListFrame.h"

ClassImp(AliFileListFrame)

//___________________________________________________________________________
AliFileListFrame::AliFileListFrame(const TGWindow *main, UInt_t w, UInt_t h): TGCompositeFrame(main, w, h), fTableLayout(0) {
  // Creates a composite frame containing a filelist widget.
   
  // use hierarchical cleaning
  SetCleanup(kDeepCleanup);
  
  fHFrame1 = new TGHorizontalFrame(this, 300, 100, kFixedWidth);
  AddFrame(fHFrame1, new TGLayoutHints(kLHintsTop|kLHintsExpandX));
  
  
  BuildQueryPathFrame();
  
  // creates the Tags for Query display
  
  fTags = new TObjArray();
  
  fTags->Add(new TObjString("type"));
  fTags->Add(new TObjString("owner"));
  fTags->Add(new TObjString("gowner"));
  fTags->Add(new TObjString("perm"));
  fTags->Add(new TObjString("size"));
  fTags->Add(new TObjString("ctime"));
  fTags->Add(new TObjString("lfn"));
  
  Pixel_t white;
  gClient->GetColorByName("white", white);
  
  fCanvas = new TGCanvas(this, 300, 300, kFixedWidth);
  fContents = new TGCompositeFrame(fCanvas->GetViewPort(), 300, 300, 
				   kSunkenFrame | kFixedWidth, white);
  
  fCanvas->SetContainer(fContents);
  AddFrame(fCanvas, new TGLayoutHints(kLHintsBottom|kLHintsExpandX, 5,5,5,0));
  
  MapSubwindows();
  MapWindow();
  Resize();
}

//___________________________________________________________________________
AliFileListFrame::~AliFileListFrame() {
  // AliFileListFrame Destructor
  // Cleanup.
  
  delete fContents;
  DeleteWindow();  // deletes fMain
}

//___________________________________________________________________________
void AliFileListFrame::SetQueryPath(const char* path) {
  // Set the Query Path
  fTextQueryPath->SetText(path);
}

//___________________________________________________________________________
const char* AliFileListFrame::GetQueryPath() {
  // Get the Query Path
  return fTextQueryPath->GetText();
}

//___________________________________________________________________________
const char* AliFileListFrame::GetQueryPattern() {
  // Get the Query Patttern
  return fTextQueryPattern->GetText();
}

//___________________________________________________________________________
void AliFileListFrame::BuildQueryPathFrame() {
  // Build the Query Path Frame
  
  fHFrame1 = new TGHorizontalFrame(this, 100, 100, kRaisedFrame);
  AddFrame(fHFrame1, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fVFrame1 = new TGVerticalFrame(fHFrame1, 100, 50);
  fHFrame1->AddFrame(fVFrame1, new TGLayoutHints(kLHintsLeft, 5,5,5,5));  
  
  fVFrame2 = new TGVerticalFrame(fHFrame1, 100, 50);
  fHFrame1->AddFrame(fVFrame2, new TGLayoutHints(kLHintsRight, 5,5,5,5));
  
  // fVFrame1, for the labels   
  
  fLabel1 = new TGLabel(fVFrame1, new TGString("Query Path"));
  fVFrame1->AddFrame(fLabel1, new TGLayoutHints(kLHintsTop, 5,5,10,5));
  
  fLabel2 = new TGLabel(fVFrame1, new TGString("Query Pattern"));
  fVFrame1->AddFrame(fLabel2, new TGLayoutHints(kLHintsCenterY, 5,5,10,5));
  
  fLabel3 = new TGLabel(fVFrame1, new TGString("Max. Results"));
  fVFrame1->AddFrame(fLabel3, new TGLayoutHints(kLHintsCenterY, 5,5,10,5));
   
  // fVFrame2 for the text boxes
  
  fTextQueryPath = new TGTextEntry(fVFrame2, new TGTextBuffer(50));
  fTextQueryPath->SetEnabled(false);
  fVFrame2->AddFrame(fTextQueryPath, 
		     new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fTextQueryPattern = new TGTextEntry(fVFrame2, new TGTextBuffer(20), 0);
  fVFrame2->AddFrame(fTextQueryPattern, 
		     new TGLayoutHints(kLHintsTop, 5,5,5,5));
  fTextQueryPattern->SetText("*.root");
  
  fNumMaxResults = new TGNumberEntry(fVFrame2, 100);
  fNumMaxResults->SetLimits(TGNumberFormat::kNELLimitMin, 0);
  fVFrame2->AddFrame(fNumMaxResults, 
		     new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fButtonRun =  new TGTextButton(fVFrame2, "          Run         ", 0);
  fVFrame2->AddFrame(fButtonRun, 
		     new TGLayoutHints(kLHintsRight, 0,5,0,0));
  
  fButtonRun->Connect("Clicked()", "AliFileListFrame", this, "RunQuery()");
}

//___________________________________________________________________________
void AliFileListFrame::RunQuery() {
  // Run Query
  
  if(!gGrid)
    return;
  
  SetCleanup(kDeepCleanup);   
  
  fContents->Cleanup();
  
  TGridResult *result = gGrid->Query(GetQueryPath(), GetQueryPattern());
  
  if(!result)
    return;

  int nrows = result->GetEntries();
  int ncols = fTags->GetEntries();
  
  if(nrows == 0){ // there is nothing to show
    TString msg = TString("The query of \"");
    (msg+=GetQueryPattern())+="\" pattern is empty ";
    new TGMsgBox(gClient->GetRoot(), this, "Empty query", msg.Data(), 0, kMBOk);
    return;
  }
  
  //if(fTableLayout != NULL) delete fTableLayout;
  fTableLayout = new TGTableLayout(fContents, nrows+1, ncols);
  fContents->SetLayoutManager(fTableLayout);
  
  Pixel_t white;
  gClient->GetColorByName("white", white);
  
  TObjString * tag;
  
  for (int i=0; i < ncols; i++) {
    for (int j=0; j < nrows; j++) {
      tag = (TObjString*) (fTags->At(i));
      TGLabel * label = 
	new TGLabel(fContents, new TGString(result->GetKey(j,tag->GetName())));        
      label->SetBackgroundColor(white);
      TGTableLayoutHints *thint = 
	new TGTableLayoutHints(i,i+1,j+1,j+2, kLHintsExpandX | kLHintsLeft, 5,30,5,5);
      fContents->AddFrame(label,thint);	
    }//for_j
  }//for_i
  
  for (int j=0; j < ncols; j++) {
    tag = (TObjString*) (fTags->At(j));
    TGTextButton * button = new TGTextButton(fContents, tag->GetName(), 0);
    TGTableLayoutHints *thint = 
      new TGTableLayoutHints(j,j+1,0,1, kLHintsFillX, 0,0,0,0);
    fContents->AddFrame(button,thint);
  }
  
  MapSubwindows();
  MapWindow();
  
  Resize();
}

//___________________________________________________________________________
void AliFileListFrame::DisplayObject(const TString& fname,const TString& name) const {
  // Browse object located in file.
  
  TDirectory *sav = gDirectory;
  
  static TFile *file = 0;
  if (file) delete file;     // close
  file = new TFile(fname);   // reopen
  
  TObject* obj = file->Get(name);
  if (obj) {
    if (!obj->IsFolder()) {
      obj->Browse(0);
    } else obj->Print();
  }
  gDirectory = sav;
}

//___________________________________________________________________________
void AliFileListFrame::OnDoubleClick(TGLVEntry *f, Int_t btn) {
  // Handle double click.
  
  if (btn != kButton1) return;
  
  // set kWatch cursor
  ULong_t cur = gVirtualX->CreateCursor(kWatch);
  gVirtualX->SetCursor(fContents->GetId(), cur);
  
  TString name(f->GetTitle());
  const char* fname = (const char*)f->GetUserData();
  
  if (fname) {
    DisplayObject(fname, name);
  } else if (name.EndsWith(".root")) {
    //      DisplayFile(name);
  } else {
    // DisplayDirectory(name);
  }
  
  // set kPointer cursor
  cur = gVirtualX->CreateCursor(kPointer);
  gVirtualX->SetCursor(fContents->GetId(), cur);
}

