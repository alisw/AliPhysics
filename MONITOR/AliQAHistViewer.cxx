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

////////////////////////////////////////////////////////////////////////////
//
//  QA histogram viewer
//
//  origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
///////////////////////////////////////////////////////////////////////////

#include "AliQAHistViewer.h"

ClassImp(AliQAHistViewer)

//_________________________________________________________________________
void AliQAHistViewer::DoDrawNext()
{
   Int_t rows = 2;
   Int_t cols = 2;
   TString oldDirStr;
   TString newDirStr;
   oldDirStr = fQANavigator->GetDirName();

   UpdateAllPathComboBoxes();

   TCanvas *c1 = fEcan->GetCanvas();
   c1->Clear();
   c1->Divide(rows,cols);
   for (Int_t i=1; i<=rows*cols;i++)
   {
       newDirStr = fQANavigator->GetDirName();
       if (newDirStr!=oldDirStr)
       {
           oldDirStr=newDirStr;
           break;
       }
       c1->cd(i);
       TH1* hist;
       if (fQANavigator->GetHistogram(hist))
       {
          if (hist) hist->Draw();
       }
       if (!fQANavigator->Next())
       {
           break;
       }
   }
   c1->Update();
}

//_________________________________________________________________________
void AliQAHistViewer::DoDrawPrev()
{
   Int_t rows = 2;
   Int_t cols = 2;
   TString oldDirStr;
   TString newDirStr;
   oldDirStr = fQANavigator->GetDirName();

   UpdateAllPathComboBoxes();

   TCanvas *c1 = fEcan->GetCanvas();
   c1->Clear();
   c1->Divide(rows,cols);
   for (Int_t i=1; i<=rows*cols;i++)
   {
       newDirStr = fQANavigator->GetDirName();
       if (newDirStr!=oldDirStr)
       {
           oldDirStr=newDirStr;
           break;
       }
       c1->cd(i);
       TH1* hist;
       if (fQANavigator->GetHistogram(hist))
       {
          if (hist) hist->Draw();
       }
       if (!fQANavigator->Prev())
       {
           break;
       }
   }
   c1->Update();
}

//_________________________________________________________________________
void AliQAHistViewer::DoExit()
{
   printf("Exit application...");
   gApplication->Terminate(0);
}

//_________________________________________________________________________
AliQAHistViewer::AliQAHistViewer(const TGWindow *p, UInt_t w, UInt_t h, Bool_t embed) :
    fEcan(NULL),
    fQANavigator(new AliQAHistNavigator()),
    TGMainFrame(p, w, h),
    fFileListBox(NULL),
    fDetectorListBox(NULL),
    fLevelListBox(NULL),
    fHistListBox(NULL),
    fExpertMode(NULL),
    fIsEmbedded(embed)
{
   //initialize the QA navigator
   // horizontal frame with comboboxes for navigation
   TGHorizontalFrame *hframenav = new TGHorizontalFrame(this, 200,40);
   fFileListBox = new TGComboBox(hframenav); 
   fFileListBox->Connect("Selected(Int_t)", "AliQAHistViewer", this, "DoSetFile(Int_t)");
   fFileListBox->Resize(150,20);
   hframenav->AddFrame(fFileListBox, new TGLayoutHints(kLHintsExpandY|kLHintsLeft));
   fDetectorListBox = new TGComboBox(hframenav);
   fDetectorListBox->Connect("Selected(Int_t)", "AliQAHistViewer", this, "DoSetDetector(Int_t)");
   fDetectorListBox->Resize(100,20);
   hframenav->AddFrame(fDetectorListBox, new TGLayoutHints(kLHintsLeft));
   fLevelListBox = new TGComboBox(hframenav);
   fLevelListBox->Connect("Selected(Int_t)", "AliQAHistViewer", this, "DoSetLevel(Int_t)");
   fLevelListBox->Resize(100,20);
   hframenav->AddFrame(fLevelListBox, new TGLayoutHints(kLHintsLeft));
   fHistListBox = new TGComboBox(hframenav);
   fHistListBox->Connect("Selected(Int_t)", "AliQAHistViewer", this, "DoSetHistogram(Int_t)");
   fHistListBox->Resize(250,20);
   hframenav->AddFrame(fHistListBox, new TGLayoutHints(kLHintsLeft));
   AddFrame(hframenav, new TGLayoutHints((kLHintsLeft|kLHintsExpandX), 5,5,5,5));
   UpdateAllPathComboBoxes();
   fExpertMode = new TGCheckButton(hframenav,"Expert");
   hframenav->AddFrame(fExpertMode,new TGLayoutHints(kLHintsLeft, 0, 4, 3, 0));
   fExpertMode->SetToolTipText("Show expert histograms");
   fExpertMode->Connect("Toggled(Bool_t)", "AliQAHistViewer", this, "DoSetExpertMode(Bool_t)");
   // Create the embedded canvas
   fEcan = new TRootEmbeddedCanvas(0,this,800,600);
   Int_t wid = fEcan->GetCanvasWindowId();
   TCanvas *myc = new TCanvas("MyCanvas", 10,10,wid);
   fEcan->AdoptCanvas(myc);
   //myc->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliQAHistViewer",this, 
   //            "EventInfo(Int_t,Int_t,Int_t,TObject*)");

   AddFrame(fEcan, new TGLayoutHints(kLHintsTop | kLHintsLeft | 
                                     kLHintsExpandX  | kLHintsExpandY,0,0,1,1));
  
   // Create a horizontal frame containing the buttons
   TGHorizontalFrame *hframebuttons = new TGHorizontalFrame(this, 200, 40); 
   TGTextButton *prev = new TGTextButton(hframebuttons, "&Prev");
   prev->Connect("Clicked()", "AliQAHistViewer", this, "DoDrawPrev()");
   hframebuttons->AddFrame(prev, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   TGTextButton *next = new TGTextButton(hframebuttons, "&Next");
   next->Connect("Clicked()", "AliQAHistViewer", this, "DoDrawNext()");
   hframebuttons->AddFrame(next, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   AddFrame(hframebuttons, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

   if ((!fIsEmbedded))
   {
       TGTextButton *exit = new TGTextButton(hframebuttons, "&Exit ");
       exit->Connect("Pressed()", "AliQAHistViewer", this, "DoExit()");
       hframebuttons->AddFrame(exit, new TGLayoutHints(kLHintsRight, 5, 5, 3, 4));
   }
   
   // Set a name to the main frame   
   SetWindowName("Quality Assurance Monitoring");
   MapSubwindows();

   // Initialize the layout algorithm via Resize()
   Resize(GetDefaultSize());

   // Map main frame
   MapWindow();
   DoDrawNext();
}

//_________________________________________________________________________
AliQAHistViewer::~AliQAHistViewer()
{
   // Clean up main frame...
   Cleanup();
   delete fEcan;
   delete fQANavigator;
}

//_________________________________________________________________________
void AliQAHistViewer::FillComboBoxWithListEntries( TGComboBox* box, const TList* list )
{
    box->RemoveAll();
    Int_t i=0;
    TIter listiter(list);
    TObject* o = NULL;
    while ((o = (TObject*)listiter.Next()))
    {
        TString name = o->GetName();
        box->AddEntry( name.Data(), i++ );
    }
}

//_________________________________________________________________________
void AliQAHistViewer::UpdateAllPathComboBoxes()
{
    if (!fQANavigator->InitOK()) return;
    FillComboBoxWithListEntries( fFileListBox, (TList*)fQANavigator->GetFileList()->GetDirs() );
    FillComboBoxWithListEntries( fDetectorListBox, (TList*)fQANavigator->GetDetectorList()->GetDirs() );
    FillComboBoxWithListEntries( fLevelListBox, (TList*)fQANavigator->GetLevelList()->GetDirs() );
    FillComboBoxWithListEntries( fHistListBox, (TList*)fQANavigator->GetItemList() );
    fFileListBox->Select(fQANavigator->GetCurrListOfFiles()->GetDirs()->IndexOf(fQANavigator->GetCurrFile()),kFALSE);
    fDetectorListBox->Select(fQANavigator->GetCurrFile()->GetDirs()->IndexOf(fQANavigator->GetCurrDetector()),kFALSE);
    fLevelListBox->Select(fQANavigator->GetCurrDetector()->GetDirs()->IndexOf(fQANavigator->GetCurrLevel()),kFALSE);
    fHistListBox->Select(fQANavigator->GetItemList()->IndexOf(fQANavigator->GetCurrItem()),kFALSE);
}

//_________________________________________________________________________
void AliQAHistViewer::DoSetFile( Int_t s )
{
    fQANavigator->SetFile(s);
    DoDrawNext();
}

//_________________________________________________________________________
void AliQAHistViewer::DoSetDetector( Int_t s )
{
    fQANavigator->SetDetector(s);
    DoDrawNext();
}

//_________________________________________________________________________
void AliQAHistViewer::DoSetLevel( Int_t s )
{
    fQANavigator->SetLevel(s);
    DoDrawNext();
}

//_________________________________________________________________________
void AliQAHistViewer::DoSetHistogram( Int_t s )
{
    fQANavigator->SetItem(s);
    DoDrawNext();
}

//_________________________________________________________________________
void AliQAHistViewer::DoSetExpertMode(Bool_t mode)
{
    fQANavigator->SetExpertMode(mode);
    DoDrawNext();
}
