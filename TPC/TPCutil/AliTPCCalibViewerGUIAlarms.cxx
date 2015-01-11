
#include <TGClient.h>
#include <TGSplitter.h>
#include <TGListTree.h>
#include <TGTextEntry.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>

#include "AliLog.h"
#include "AliTPCCalibQAChecker.h"
#include "AliTPCCalibViewerGUI.h"
#include "AliTPCCalibViewerGUItime.h"

#include "AliTPCCalibViewerGUIAlarms.h"


ClassImp(AliTPCCalibViewerGUIAlarms)

AliTPCCalibViewerGUIAlarms::AliTPCCalibViewerGUIAlarms(const TGWindow *p, UInt_t w, UInt_t h) :
  TGCompositeFrame(p,w,h),
  fCalibChecker(0x0),
  fAlarmTree(0x0),
  fMainCanvas(0x0),
  fTreeCanvas(0x0),
  fAlarmText(0x0),
  fCalibViewerGUI(0x0),
  fCalibViewerGUItime(0x0)
{
  //
  //
  //
  DrawGUI(p,w,h);
}
//______________________________________________________________________________
AliTPCCalibViewerGUIAlarms::~AliTPCCalibViewerGUIAlarms()
{
  //
  //
  //
//   gClient->FreePicture("interrupt.xpm");
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::DrawGUI(const TGWindow */*p*/, UInt_t w, UInt_t h)
{
  //
  // draw the GUI
  //

  //GUI elements
  SetCleanup(kDeepCleanup);
  
  // *****************************************************************************
  // ************************* content of this MainFrame *************************
  // *****************************************************************************
  // top level container with horizontal layout
  
  TGCompositeFrame *contLCR = new TGCompositeFrame(this, w, h, kHorizontalFrame | kFixedWidth | kFixedHeight);
  AddFrame(contLCR, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  
  // ***********************************************************************
  // ************************* content of contLCR *************************
  // ***********************************************************************
  // left container
  TGCompositeFrame *contLeft = new TGCompositeFrame(contLCR, 200, 200, kVerticalFrame | kFixedWidth | kFitHeight);
  contLCR->AddFrame(contLeft, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY, 5, 3, 3, 3));
  
   // left vertical splitter
  TGVSplitter *splitLeft = new TGVSplitter(contLCR);
  splitLeft->SetFrame(contLeft, kTRUE);
  contLCR->AddFrame(splitLeft, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0, 0, 0, 0));
  
   // right container
  TGCompositeFrame *contRight = new TGCompositeFrame(contLCR, 150, 200, kVerticalFrame | kFixedWidth | kFitHeight);
  contLCR->AddFrame(contRight, new TGLayoutHints(kLHintsTop | kLHintsRight | kLHintsExpandY, 3, 5, 3, 3));
  
   // center container
  TGCompositeFrame *contCenter = new TGCompositeFrame(contLCR, 200, 200, kVerticalFrame | kFixedWidth | kFitHeight);
  contLCR->AddFrame(contCenter, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  
   // right vertical splitter
  TGVSplitter *splitRight = new TGVSplitter(contLCR);
  splitRight->SetFrame(contRight, kFALSE);
  contLCR->AddFrame(splitRight, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0, 0, 0, 0));

  // ***********************************************************************
  // *******************     left container         ************************
  // ***********************************************************************

  //Alarm list view
  TGCanvas *treeCanvas = new TGCanvas(contLeft, 200, 200);
  fAlarmTree = new TGListTree(treeCanvas, kHorizontalFrame);
  contLeft->AddFrame(treeCanvas,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY | kLHintsTop));
  fAlarmTree->Connect("DoubleClicked(TGListTreeItem*,Int_t)","AliTPCCalibViewerGUIAlarms",this,
                     "OnDoubleClick(TGListTreeItem*,Int_t)");
  fAlarmTree->Connect("Clicked(TGListTreeItem*,Int_t)","AliTPCCalibViewerGUIAlarms",this,
                     "OnClick(TGListTreeItem*,Int_t)");
  fAlarmTree->SetColorMode(TGListTree::kColorBox);

  // ***********************************************************************
  // *******************       center container     ************************
  // ***********************************************************************

  // top container
  TGCompositeFrame *contCenterTop = new TGCompositeFrame(contCenter, 200, 200, kVerticalFrame | kFixedWidth | kFitHeight);
  contCenter->AddFrame(contCenterTop, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  
  // bottom container
  TGCompositeFrame *contCenterBottom = new TGCompositeFrame(contCenter, 200, 100, kVerticalFrame | kFixedWidth | kFixedHeight);
  contCenter->AddFrame(contCenterBottom, new TGLayoutHints( kLHintsExpandX, 0, 0, 0, 0));
  
  //--------------
  // canvas on top
  TRootEmbeddedCanvas *cEmbed = new TRootEmbeddedCanvas("Alarm_Canvas", contCenterTop, 200, 200, kFitWidth | kFitHeight);
  contCenterTop->AddFrame(cEmbed, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
//   cEmbed->GetCanvas()->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)", "AliTPCCalibViewerGUIAlarms", this, "MouseMove(Int_t, Int_t, Int_t, TObject*)");
  cEmbed->GetCanvas()->SetToolTipText("Alarm histograms are displayed in this region.");
  fMainCanvas=cEmbed->GetCanvas();
  //--------------
  // alarm text on bottom
  // canvas
  TGCanvas *textCanvas = new TGCanvas(contCenterBottom, 200, 200);
  fAlarmText = new TGLabel(textCanvas->GetViewPort(),"Alarm descriptions can be found here");
  textCanvas->SetContainer(fAlarmText);
  contCenterBottom->AddFrame(textCanvas,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY | kLHintsTop));
  
//   textCanvas->AddFrame(fAlarmText,new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 0, 0));
  
  
  // ***********************************************************************
  // *******************       right container      ************************
  // ***********************************************************************
  
  TGGroupFrame *grpInfo = new TGGroupFrame(contRight, "Alarm Info", kVerticalFrame | kFitWidth | kFitHeight);
  contRight->AddFrame(grpInfo, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
  
  
  SetWindowName("AliTPCCalibViewer GUI - Alarms");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
  
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::AddSubItems(AliTPCCalibQAChecker *checker, TGListTreeItem *item)
{
  //
  //
  //
  TGListTreeItem *newitem=fAlarmTree->AddItem(item,checker->GetName(),(void*)checker);
  newitem->SetColor(checker->GetQualityColor());
  newitem->SetTipText(checker->GetTitle());
  if (checker->HasSubCheckers()) {
    AliTPCCalibQAChecker *ch=0x0;
    while ( (ch=checker->NextSubChecker()) ) AddSubItems(ch,newitem);
  } else {
    newitem->SetPictures(gClient->GetPicture("interrupt.xpm"),gClient->GetPicture("interrupt.xpm"));
  }
  
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::InitBrowser()
{
  //
  //
  //
  if (!fAlarmTree){
    AliError("Alarms not set!");
    return;
  }
  AddSubItems(fCalibChecker,0);
  fAlarmTree->ClearViewPort();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::UpdateSubItem(TGListTreeItem *item)
{
  //
  //
  //

//   printf("name: %s\n", item->GetText());
  AliTPCCalibQAChecker *checker=dynamic_cast<AliTPCCalibQAChecker*>((TObject*)item->GetUserData());
  if (checker){
    item->SetColor(checker->GetQualityColor());
  } else {
    item->ClearColor();
  }
  TGListTreeItem *nextItem=0x0;
  if ( (nextItem=item->GetFirstChild())  ) UpdateSubItem(nextItem);
  if ( (nextItem=item->GetNextSibling()) ) UpdateSubItem(nextItem);
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::UpdateBrowser()
{
  //
  //
  //
  UpdateSubItem(fAlarmTree->GetFirstItem());
  fAlarmTree->ClearViewPort();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::ResetBrowser()
{
  //
  //
  //
  if (!fAlarmTree->GetFirstItem()) return;
  fAlarmTree->DeleteChildren(fAlarmTree->GetFirstItem());
  fAlarmTree->DeleteItem(fAlarmTree->GetFirstItem());
  fAlarmTree->ClearViewPort();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::OpenSubItems(TGListTreeItem *item)
{
  //
  //
  //
  fAlarmTree->OpenItem(item);
  TGListTreeItem *nextItem=0x0;
  if ( (nextItem=item->GetFirstChild())  ) OpenSubItems(nextItem);
  if ( (nextItem=item->GetNextSibling()) ) OpenSubItems(nextItem);
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::OpenAllItems()
{
  //
  //
  //
  if (!fAlarmTree->GetFirstItem()) return;
  OpenSubItems(fAlarmTree->GetFirstItem());
  fAlarmTree->ClearViewPort();
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::OnDoubleClick(TGListTreeItem* item, Int_t /*id*/)
{
  //
  //
  //
//   printf("DoubleClicked");
  if (!fCalibViewerGUI && !fCalibViewerGUItime) return;
  TGTextEntry *draw=0x0;
  TGTextEntry *cuts=0x0;
  TGTextEntry *opt=0x0;
  if (fCalibViewerGUI){
    draw=fCalibViewerGUI->GetDrawEntry();
    cuts=fCalibViewerGUI->GetCutsEntry();
    opt=fCalibViewerGUI->GetDrawOptEntry();
  } else if (fCalibViewerGUItime){
    draw=fCalibViewerGUItime->GetDrawEntry();
    cuts=fCalibViewerGUItime->GetCutsEntry();
    opt=fCalibViewerGUItime->GetDrawOptEntry();
  }
  draw->SetText(((AliTPCCalibQAChecker*)item->GetUserData())->GetDrawString());
  cuts->SetText(((AliTPCCalibQAChecker*)item->GetUserData())->GetCutsString());
  opt->SetText(((AliTPCCalibQAChecker*)item->GetUserData())->GetDrawOptString());
}
//______________________________________________________________________________
void AliTPCCalibViewerGUIAlarms::OnClick(TGListTreeItem* item, Int_t /*id*/)
{
  //
  //
  //
//   printf("Clicked");
  
  TVirtualPad *pad=gPad;
  fMainCanvas->cd();
  fMainCanvas->Clear();
  AliTPCCalibQAChecker *checker=(AliTPCCalibQAChecker*)item->GetUserData();
  checker->Draw("nobc");
  TString text(checker->GetQualityDescription());
  if (text.IsNull()) fAlarmText->SetText("No description available!");
  else fAlarmText->SetText(text.Data());
  gPad->Modified();
  gPad->Update();
  if (pad) pad->cd();
  
}
//______________________________________________________________________________
AliTPCCalibViewerGUIAlarms* AliTPCCalibViewerGUIAlarms::Show()
{
  //
  //
  //
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("AliTPCCalibViewer GUI");
  frmMain->SetCleanup(kDeepCleanup);

  AliTPCCalibViewerGUIAlarms* calibViewer1 = new AliTPCCalibViewerGUIAlarms(frmMain, 1000, 600);
  frmMain->AddFrame(calibViewer1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();

  return calibViewer1;
}
//______________________________________________________________________________

