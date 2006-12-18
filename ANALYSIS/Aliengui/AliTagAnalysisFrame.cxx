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
//           AliTagAnalysisFrame class
//   The class that deals with the event tag tab of the GUI
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TGMsgBox.h"
#include "TGListBox.h"
#include "TGComboBox.h"

#include "TSystem.h"
#include "TChain.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TEventList.h"

#include "AliRunTagCuts.h"
#include "AliEventTagCuts.h"
#include "AliTagAnalysis.h"

#include "AliAnalysisGUI.h"
#include "AliAlienBrowser.h"
#include "AliTagFrame.h"
#include "AliTagAnalysisFrame.h"

ClassImp(AliTagAnalysisFrame)

//___________________________________________________________________________
AliTagAnalysisFrame::AliTagAnalysisFrame(const TGWindow *main, UInt_t w, UInt_t h, AliAnalysisGUI* fAliAnalysisGUI): 
  TGMainFrame(main, w, h, kHorizontalFrame), 
  fkNumberOfTags(3), 
  fVFrame1(0), fVFrame2(0),
  fGroup1(0), fGroup2(0), fGroup3(0),
  fAliAnalysisGUI(fAliAnalysisGUI),
  fTagFrame(0), fAliEnBrowser(0),
  fLocalLabel1(0), fLocalPath(0),
  fLocalButton(0), fButtonInsert(0), fButtonRun(0),
  fComboEventTagCut(0), fGridLabel1(0),
  fGridPath(0), fGridButton(0), fButtonInsert2(0), fButtonRun2(0),
  fComboEventTagCut2(0), fTagResult(0),
  fAnalysisChain(0), fListBox(0),
  fBrowser(NULL), fBrowserButton(NULL),
  fAliTagAnalysis(0), fAliRunCuts(0),
  fAliEventCuts(0), fEventTagCutsName(0) {
   // Constructor.

/*
   // lazy initialization to fEventTagCutsName
   const char *tmp[] ={ "Vx", "Vy", "Vz", "Participants", "Impact parameter", "Primary vertex",
		       "ZDC - neutron 1", "ZDC - proton 1", "ZDC - neutron 2", "ZDC - proton 2",
		       "ZDC EM", "TO VertexZ",
		       "Multiplicity", "Positive Multiplicity", "Negative Multiplicity", 
		       "Neutral Multiplicity", "VO", "Cascades", "Kinks", 
		       "Jet Energy", "Hard Photons Candidates", "Neutral Energy", 
		       "Charged above 1 GeV", "Charged above 3 GeV", "Charged above 10 GeV",
		       "Muons above 1 GeV", "Muons above 3 GeV", "Muons above 10 GeV", 
		       "Electron above 1 GeV", "Electron above 3 GeV", "Electron above 10 GeV",
		       "Electrons range", "Muons range", "Pions range", "Kaons range", 
		       "Protons range", "Lambda range", "Photons range", "PiOs range", 
		       "Neutrons range", "KaonOs range"
   };
 */

  const char *tmp[] = {"MultiplicityRange","VOsRange", "NPionRange" };
  fEventTagCutsName = tmp;
  
  //   fEventTagCutsName = new TList();
  
  // fEventTagCutsName[0] = "NegMultiplicityRange";
  //    fEventTagCutsName[1] = "VOsRange";
  //   fEventTagCutsName[2] = "NPionRange";

  fVFrame1 = new TGVerticalFrame(this, 200, 150);
  this->AddFrame(fVFrame1, new TGLayoutHints(kLHintsLeft, 5,5,5,5));
  
  //  Local Group
  fGroup1 = new TGGroupFrame(fVFrame1, "Local", kVerticalFrame);
  fGroup1->SetTitlePos(TGGroupFrame::kLeft); // left aligned
  fVFrame1->AddFrame(fGroup1, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  BuildLocalGroup(fGroup1);
  
  //  Grid Group
  fGroup2 = new TGGroupFrame(fVFrame1, "Grid", kVerticalFrame);
  fGroup2->SetTitlePos(TGGroupFrame::kLeft); // left aligned
  fVFrame1->AddFrame(fGroup2, new TGLayoutHints(kLHintsBottom, 5,5,5,5));
  
  BuildGridGroup(fGroup2);
  
  // Vertical Frame 2
  
  fVFrame2 = new TGVerticalFrame(this, 200, 200);
  AddFrame(fVFrame2, new TGLayoutHints(kLHintsRight| kLHintsExpandX 
				       | kLHintsExpandY, 5,5,5,5));
  
  fGroup3 = new TGGroupFrame(fVFrame2, "Results", 
			     kVerticalFrame | kFitWidth | kFitHeight);
  fGroup3->SetTitlePos(TGGroupFrame::kLeft); // left aligned
  fVFrame2->AddFrame(fGroup3, 
		     new TGLayoutHints(kLHintsTop | kLHintsExpandX 
				       | kLHintsExpandY, 5,5,5,5));      
  
  fListBox = new TGListBox(fGroup3); 
  fGroup3->AddFrame(fListBox, 
		    new TGLayoutHints(kLHintsTop | kLHintsExpandX |
				      kLHintsExpandY, 5,5,5,5));  
  
  fAliTagAnalysis = new AliTagAnalysis(); 
  fAliRunCuts = new AliRunTagCuts();
  fAliEventCuts = new AliEventTagCuts();
  
  MapSubwindows();
  Resize();
  MapWindow();
}

//___________________________________________________________________________
AliTagAnalysisFrame::~AliTagAnalysisFrame() {
  // AliTagAnalysisFrame dctor.
  
  delete fGroup1;
  delete fLocalLabel1;
  delete fLocalPath;
  delete fLocalButton;
  delete fGroup2;
  delete fGridLabel1;
  delete fGridPath;
  delete fGridButton;
  
  delete fAliTagAnalysis;
  delete fAliRunCuts;
  delete fAliEventCuts;
  delete fTagResult;
  delete fAnalysisChain;
  
  delete fTagFrame;
}

//___________________________________________________________________________
void AliTagAnalysisFrame::AddResult (const char* line) {
  // Add a new line in the result group box.
  
  //    fGroup3->AddFrame(new TGLabel(fGroup3, new TGString(line)), 
  // 		     new TGLayoutHints(kLHintsTop, 5,5,5,5));  
  
  fListBox->AddEntry(line, fListBox->GetNumberOfEntries()); 
  
  MapSubwindows();
  Resize();
  MapWindow();   
}

//___________________________________________________________________________
void AliTagAnalysisFrame::BuildLocalGroup (TGCompositeFrame* frame) {
  // The Local Group Frame
  fLocalLabel1 = new TGLabel(frame, new TGString("Chain Local Tag Path"));
  frame->AddFrame(fLocalLabel1, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fLocalPath = new TGTextEntry(frame, new TGTextBuffer(40));
  fLocalPath->SetEnabled(false);
  frame->AddFrame(fLocalPath, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fLocalButton = new TGTextButton(frame, "Browse...", 0);
  frame->AddFrame(fLocalButton, new TGLayoutHints(kLHintsLeft, 5,5,5,5));  
 
  fLocalButton->Connect("Clicked()", "AliTagAnalysisFrame", this, "LocalBrowse()");
 
  fComboEventTagCut = new TGComboBox(frame, "Select Tag Cuts...", 1);
  frame->AddFrame(fComboEventTagCut,
		  new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5,5,5,5));
  
  for(int i=0; i!=fkNumberOfTags; i++)
    fComboEventTagCut->AddEntry(fEventTagCutsName[i],i);
  
  fComboEventTagCut->Resize(150, 20);
  
  fButtonInsert = new TGTextButton(frame, "Insert Tag Cuts Range", 2);
  frame->AddFrame(fButtonInsert,
		  new TGLayoutHints(kLHintsLeft | kLHintsTop, 5,5,5,5));
  
  fButtonInsert->Connect("Clicked()", "AliTagAnalysisFrame", this,
			 "InsertTagCutsRangeLocal()");
  
  fButtonRun = new TGTextButton(frame,     "         Run        ", 3);
  frame->AddFrame(fButtonRun,
		  new TGLayoutHints(kLHintsTop | kLHintsRight, 5,5,5,5));
  
  fButtonRun->Connect("Clicked()", "AliTagAnalysisFrame", this, "RunLocal()"); 
}

//___________________________________________________________________________
void AliTagAnalysisFrame::BuildGridGroup (TGCompositeFrame* frame) {
  // The Grid Group Frame
  
  fGridLabel1 = new TGLabel(frame, new TGString("Chain Grid Tag Path"));
  frame->AddFrame(fGridLabel1, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fGridPath = new TGTextEntry(frame, new TGTextBuffer(40));
  fGridPath->SetEnabled(false);
  //   fGridPath->SetText("/alice/cern.ch/user/p/pchrista/PDC06/Tags/pp/1");
  frame->AddFrame(fGridPath, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fGridButton = new TGTextButton(frame, "Browse...", 0);
  frame->AddFrame(fGridButton, new TGLayoutHints(kLHintsLeft, 5,5,5,5));
  
  fGridButton->Connect("Clicked()", "AliTagAnalysisFrame", this, "GridBrowse()");
  
  fComboEventTagCut2 = new TGComboBox(frame, "Select Tag Cuts...", 1);
  frame->AddFrame(fComboEventTagCut2, 
		      new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5,5,5,5));
  
  for(int i=0; i!=fkNumberOfTags; i++)
    fComboEventTagCut2->AddEntry(fEventTagCutsName[i],i);
  
  fComboEventTagCut2->Resize(150, 20);
  
  fButtonInsert2 = new TGTextButton(frame, "Insert Tag Cuts Range", 2);
  frame->AddFrame(fButtonInsert2, 
		  new TGLayoutHints(kLHintsLeft, 5,5,5,5));
  
  fButtonInsert2->Connect("Clicked()", "AliTagAnalysisFrame", this, 
			  "InsertTagCutsRangeGrid()");
  
  fButtonRun2 = new TGTextButton(frame,"         Run        " , 2);
  frame->AddFrame(fButtonRun2, 
		  new TGLayoutHints(kLHintsRight, 5,5,5,5));
  
  fButtonRun2->Connect("Clicked()", "AliTagAnalysisFrame", this, "RunGrid()");
  
}

//___________________________________________________________________________
void AliTagAnalysisFrame::LocalBrowse() {
  // Browse local directories.
  
  fBrowser = new TGTransientFrame(gClient->GetRoot(), fAliAnalysisGUI, 450, 200);
  fAliEnBrowser = new AliAlienBrowser(fBrowser, 300, 200, this, 
				      "AliTagAnalysisFrame", kLocalBrowse);
  fBrowser->AddFrame(fAliEnBrowser, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  fBrowserButton = new TGTextButton(fBrowser, "  OK  ", 0);
  fBrowser->AddFrame(fBrowserButton, new TGLayoutHints(kLHintsRight, 5,5,5,5));
  fBrowserButton->Connect("Clicked()", "AliTagAnalysisFrame", this, "OnOKButton()");
  
  fAliEnBrowser->AddItem(0, "/");
  
  fAliEnBrowser->GotoDir(gSystem->pwd());
  
  fBrowser->MapSubwindows();
  fBrowser->Resize();
  fBrowser->MapWindow();
}

//___________________________________________________________________________
void AliTagAnalysisFrame::GridBrowse() {
  // Opens a browser for grid directories.
  
  if (!fAliAnalysisGUI->IsConnected()){
    new TGMsgBox(gClient->GetRoot(), this, "Connect", 
		 "Please connect to AliEn", 0, kMBOk);
    return;
  }
  
  fBrowser = new TGTransientFrame(gClient->GetRoot(), fAliAnalysisGUI, 450, 200);
  
  fAliEnBrowser = new AliAlienBrowser(fBrowser, 300, 200, this, 
				      "AliTagAnalysisFrame", kGridBrowse);
  fBrowser->AddFrame(fAliEnBrowser, new TGLayoutHints(kLHintsTop, 5,5,5,5));
  
  fBrowserButton = new TGTextButton(fBrowser, "  OK  ", 0);
  fBrowser->AddFrame(fBrowserButton, new TGLayoutHints(kLHintsRight, 5,5,5,5));
  
  fBrowserButton->Connect("Clicked()", "AliTagAnalysisFrame", this, "OnOKButton()");
  
  fAliEnBrowser->AddItem(0, "/");
  
  fAliEnBrowser->GotoDir(gGrid->GetHomeDirectory());
  
  fBrowser->MapSubwindows();
  fBrowser->Resize();
  fBrowser->MapWindow();
}

//___________________________________________________________________________
void AliTagAnalysisFrame::InsertTagCutsRangeLocal() {
  // slot
  InsertTagCutsRange(fComboEventTagCut->GetSelected());
}


//___________________________________________________________________________
void AliTagAnalysisFrame::InsertTagCutsRangeGrid() {
  // slot
  InsertTagCutsRange(fComboEventTagCut2->GetSelected());
}

//___________________________________________________________________________
void AliTagAnalysisFrame::InsertTagCutsRange(Int_t id) {
   // insert the event tag range

   // if nth is selected
  if(id == -1)
    return;
  
  
  switch(id){
  case 0: // SetMultiplicity Range
    
    fTagFrame = new AliTagFrame(gClient->GetRoot(), this, 400, 200, kHorizontalFrame, fComboEventTagCut->GetTextEntry()->GetText(), fComboEventTagCut->GetSelected(), kRangeMinMax);
    
    Int_t min = fTagFrame->GetRangeMin();
    Int_t max = fTagFrame->GetRangeMax();
    
    fAliEventCuts->SetMultiplicityRange(min, max);
    
    TString res = TString("Multiplicity Range Min: ");
    res += min;
    res += " Max: ";
    res += max;
    
    AddResult(res.Data());
    
    break;
  }
}

//___________________________________________________________________________
void AliTagAnalysisFrame::RunLocal() {
  // Run local query
#ifdef GUIDEBUG     
  printf("*******************************\n");
  printf("*** Querying the tags       ***\n");
  printf("*******************************\n");
#endif
  
  //local tags
  fAliTagAnalysis->ChainLocalTags(fLocalPath->GetText());
  
#ifdef GUIDEBUG     
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
#endif   
  
  fAnalysisChain = new TChain("esdTree");
  fAnalysisChain = fAliTagAnalysis->QueryTags(fAliRunCuts,fAliEventCuts);
  
  TString res = TString("Number of Accepted Events: ");
  res += fAnalysisChain->GetEventList()->GetN();
  
  AddResult(res.Data()); 
}

//___________________________________________________________________________
void AliTagAnalysisFrame::RunGrid() {
  // Run Grid query
  
  
  //   fGroup3->SetCleanup(kDeepCleanup);
  
  if (!fAliAnalysisGUI->IsConnected()){
    new TGMsgBox(gClient->GetRoot(), this, "Connect", 
		 "Please connect to AliEn", 0, kMBOk);
    return;
  }
  
  
#ifdef GUIDEBUG 
  printf("*******************************\n");
  printf("*** Querying the tags       ***\n");
  printf("*******************************\n");
#endif
  
  
  //   TGridResult* TagResult = gGrid->Query("/alice/cern.ch/user/p/pchrista/PDC06/Tags/pp/1","*tag.root","","");
  fTagResult = gGrid->Query(fGridPath->GetText(), "*tag.root", "", "");
  
  //   fAliTagAnalysis->ChainLocalTags("../tags");
  
  fAliTagAnalysis->ChainGridTags(fTagResult);
  
  //////////////////////////////////////////////////////////////////
  //Get the chain
#ifdef GUIDEBUG 
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
#endif
  
  fAnalysisChain = new TChain("esdTree");
  fAnalysisChain = fAliTagAnalysis->QueryTags(fAliRunCuts,fAliEventCuts);
  
  TString res = TString("Number of Accepted Events: ");
  res += fAnalysisChain->GetEventList()->GetN();
  
  AddResult(res.Data());
}

//___________________________________________________________________________
void AliTagAnalysisFrame::ProcessSelector(const char* selectorfile) {
  // Process selector
  
#ifdef GUIDEBUG      
  printf("*******************************\n");
  printf("*** Run Analysis Selector %s\n",selectorfile);
  printf("*******************************\n");
  
#endif
  
  fAnalysisChain->Process(selectorfile);
}


//___________________________________________________________________________
void AliTagAnalysisFrame::OnDoubleClick(TGListTreeItem* item, Int_t btn) {
  // Slot for double clicking.
  
  fAliEnBrowser->OnDoubleClick(item, btn);
  
  //"/alice/cern.ch/user/p/pchrista/PDC06/Tags/pp/1");
  
}

//___________________________________________________________________________
void AliTagAnalysisFrame::OnOKButton() {
  // Slot for OK button in the Transient Frame.
  
  if(fAliEnBrowser->GetBrowseType() == kLocalBrowse)
    fLocalPath->SetText(fAliEnBrowser->GetPath());  
  else if(fAliEnBrowser->GetBrowseType() == kGridBrowse)
    fGridPath->SetText(fAliEnBrowser->GetPath());  
  
  TTimer::SingleShot(150, "AliTagAnalysisFrame", fBrowser, "CloseWindow()");
}
