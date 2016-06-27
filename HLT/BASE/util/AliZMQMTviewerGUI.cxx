#include "TApplication.h"
#include "TGLabel.h"
#include "TGClient.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TPRegexp.h"
#include "TGButton.h"
#include "TGColorSelect.h"
#include "TRootEmbeddedCanvas.h"
#include "TH1F.h"
#include "TGTextEntry.h"
#include "TThread.h"
#include "AliZMQMTviewerGUI.h"
#include "AliZMQhistViewer.h"
#include "AliZMQhelpers.h"
#include "zmq.h"
#include <cstdlib>
#include <string>
#include "TTimer.h"
#include "Buttons.h"
#include "TRootCanvas.h"
#include "TGStatusBar.h"
#include "TPad.h"

ClassImpQ(AliZMQMTviewerGUI)

using namespace std;

const char* AliZMQMTviewerGUI::fUSAGE =
"ZMQhstViewer: Draw() all ROOT drawables in a message\n"
"options: \n"
" -in : data in\n"
" -sleep : how long to sleep in between requests for data in s (if applicable)\n"
" -timeout : how long to wait for the server to reply (s)\n"
" -Verbose : be verbose\n"
" -select : only show selected histograms (by regexp)\n"
" -unselect : as select, only inverted\n"
" -drawoptions : what draw option to use\n"
" -file : dump input to file and exit\n"
" -log[xyz] : use log scale on [xyz] dimension\n"
" -histstats : histogram stat box options (default 0)\n"
" -AllowResetAtSOR : 0/1 to reset at change of run\n"
;

//______________________________________________________________________
AliZMQMTviewerGUI::AliZMQMTviewerGUI(const TGWindow *p,UInt_t w,UInt_t h, int argc, char** argv)
  : TGMainFrame(p,w,h)
  , fECanvas(NULL)
  , fStatusBar(NULL)
  , fPullButton(NULL)
  , fSelectionEntry(NULL)
  , fUnSelectionEntry(NULL)
  , fThread(NULL)
  , fCanvas(NULL)
  , fThreadArgs()
  , fViewer(NULL)
  , fTimer(NULL)
  , fWindowTitle()
  , fZMQviewerConfig()
  , fInitStatus(0)
  , fSelection(NULL)
  , fUnSelection(NULL)
{
   fViewer = new AliZMQhistViewer();
   int viewerOptionRC = fViewer->ProcessOptionString(argc,argv);
   if (viewerOptionRC<1) {
     fInitStatus = -1;
     return;
   }
   fSelection = fViewer->GetSelection();
   fUnSelection = fViewer->GetUnSelection();

   // Creates widgets of the example
   SetCleanup(kDeepCleanup);
   fECanvas = new TRootEmbeddedCanvas ("Ecanvas", this, 600, 400);
   Int_t wid = fECanvas->GetCanvasWindowId();
   fCanvas = new TCanvas("", 10,10,wid);
   fECanvas->AdoptCanvas(fCanvas);   
   //fCanvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliZMQMTviewerGUI",this, 
   //            "EventInfo(Int_t,Int_t,Int_t,TObject*)");
   fCanvas->Connect("Selected(TVirtualPad*,TObject*,Int_t)","AliZMQMTviewerGUI",this, 
               "PadSelected(TVirtualPad*,TObject*,Int_t)");
   fCanvas->Connect("Picked(TPad*,TObject*,Int_t)","AliZMQMTviewerGUI",this, 
               "PadPicked(TPad*,TObject*,Int_t)");

   AddFrame(fECanvas, new TGLayoutHints(kLHintsTop | kLHintsLeft | 
                                     kLHintsExpandX  | kLHintsExpandY,0,0,1,1));
   
   
   ///////////////////////////////
   // selection / unselection frame
   TGCompositeFrame* fSelectionEntryFrame =
     new TGCompositeFrame(this, 60, 20, kHorizontalFrame | kSunkenFrame);
   
   TGLabel *selectionLabel = new TGLabel(fSelectionEntryFrame, "select:");
   fSelectionEntryFrame->AddFrame(selectionLabel,
       new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 2, 2, 2));

   fSelectionEntry = new TGTextEntry(fSelectionEntryFrame);
   fSelectionEntry->SetToolTipText("Enter a selection regexp");
   //fSelectionEntry->Resize(300, fSelectionEntry->GetDefaultHeight());
   fSelectionEntryFrame->AddFrame(fSelectionEntry,
       new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 2, 2, 2));
   fSelectionEntry->Connect("TabPressed()","AliZMQMTviewerGUI",this,"DoSelectionEntry()");
   fSelectionEntry->Connect("ReturnPressed()","AliZMQMTviewerGUI",this,"DoSelectionEntry()");
   if (fSelection) fSelectionEntry->SetText(fSelection->GetPattern().Data(), kFALSE);
   
   TGLabel *unSelectionLabel = new TGLabel(fSelectionEntryFrame, "unselect:");
   fSelectionEntryFrame->AddFrame(unSelectionLabel,
       new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 2, 2, 2));

   fUnSelectionEntry = new TGTextEntry(fSelectionEntryFrame);
   fUnSelectionEntry->SetToolTipText("this will be inversely matched");
   //fUnSelectionEntry->Resize(300, fUnSelectionEntry->GetDefaultHeight());
   fSelectionEntryFrame->AddFrame(fUnSelectionEntry,
       new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 2, 2, 2));
   fUnSelectionEntry->Connect("TabPressed()","AliZMQMTviewerGUI",this,"DoSelectionEntry()");
   fUnSelectionEntry->Connect("ReturnPressed()","AliZMQMTviewerGUI",this,"DoSelectionEntry()");
   if (fUnSelection) fUnSelectionEntry->SetText(fUnSelection->GetPattern().Data(), kFALSE);

   fPullButton = new TGTextButton(fSelectionEntryFrame, "Pull", 150);
   fPullButton->Connect("Clicked()", "AliZMQMTviewerGUI", this, "DoPullButton()");
   fPullButton->SetToolTipText("force update");
   fSelectionEntryFrame->AddFrame(fPullButton,
       new TGLayoutHints(kLHintsTop | kLHintsRight, 10, 5, 2, 2));
   AddFrame(fSelectionEntryFrame, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 1, 1, 1, 1));
   
   ///////////////////////////////
   // Create status frame containing a button and a text entry widget
   TGCompositeFrame* fStatusFrame =
     new TGCompositeFrame(this, 60, 20, kHorizontalFrame | kSunkenFrame);

   Int_t parts[] = {70, 30};
   fStatusBar = new TGStatusBar(fStatusFrame, 100, 1, kHorizontalFrame);
   fStatusBar->SetParts(parts, 2);
   fStatusBar->Draw3DCorner(kFALSE);
   fStatusFrame->AddFrame(fStatusBar,
       new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX, 5, 5, 2, 2));

   AddFrame(fStatusFrame, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 1, 1, 1, 1));

   
   int rc = alizmq_socket_init(fZMQviewerConfig, alizmq_context(), "PUSH@inproc://viewerConfig");
   if (rc < 0) {
     printf("WARNING: viewer ZMQ config socket init error: %i %s\n", rc , zmq_strerror(errno));
   }

   fViewer->SetCanvas(fCanvas);
   bool updateCanvas = false;
   fViewer->GetUpdateCanvas(&updateCanvas);

   fTimer = new TTimer(this, 1000);
   fTimer->Reset();
   fTimer->TurnOn();

   // Sets window name and shows the main frame
   fWindowTitle = "HLT hist viewer";
   SetWindowName(fWindowTitle.c_str());
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();

   // set up the threads
   fThread = new TThread("MyThread1", &AliZMQMTviewerGUI::RunViewer, (void*)this);
   fThreadArgs.fRun = kTRUE;
   fThreadArgs.fThread = fThread;
   fThreadArgs.fCanvas = fCanvas;
   fThreadArgs.viewer = fViewer;
   fThread->Run(&fThreadArgs);
}

//______________________________________________________________________________
void AliZMQMTviewerGUI::CloseWindow()
{
   // terminate and delete the thread, then quit the application
   fThreadArgs.fRun = kFALSE;
   bool term = true;
   zmq_close(fZMQviewerConfig);
   zmq_ctx_term(alizmq_context());
   fViewer->GetTerminated(&term);
   while (fThread->GetState() == TThread::kRunningState) {
      gSystem->ProcessEvents();
   }
   fThread->Join();
   delete fThread;
   TGMainFrame::CloseWindow();
   gApplication->Terminate();
}

//______________________________________________________________________________
void *AliZMQMTviewerGUI::RunViewer(void *ptr)
{
   // this method is called by the thread

   ThreadArgs_t *args = (ThreadArgs_t *)ptr;
   if (!args) return 0;
   AliZMQhistViewer* viewer = args->viewer;

   viewer->Run(NULL);

   return NULL;
}

//______________________________________________________________________________
void AliZMQMTviewerGUI::UpdateCanvas()
{
   std::string info = fViewer->GetInfo();
   if (info.size()>0) {
     info=", "+info;
   }
   std::string title = fWindowTitle + info;
   SetWindowName(title.c_str());
   fViewer->UpdateCanvas(fCanvas, fSelection, fUnSelection);
}

//______________________________________________________________________________
Bool_t AliZMQMTviewerGUI::HandleTimer(TTimer *t)
{
   UpdateCanvas();
   fStatusBar->SetText(Form("source: %s",fViewer->GetZMQconfig().c_str()), 0);
   fTimer->Reset();
   return kTRUE;
}

//______________________________________________________________________________
void AliZMQMTviewerGUI::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected)
{
  if (event==kMouseMotion) return;
  if (!selected) printf("event info null pointer\n");
  if (event==kMouseEnter) {
    fStatusBar->SetText(Form("%s",selected->GetName()),1);
  }
}

//______________________________________________________________________________
void AliZMQMTviewerGUI::PadSelected(TVirtualPad* pad, TObject *object, Int_t event)
{
  if (event==1 && pad) {
    if (pad==fCanvas) { return; }
    pad->cd();
    std::string name = "";
    std::string title = "";
    Int_t logx = pad->GetLogx();
    Int_t logy = pad->GetLogy();
    Int_t logz = pad->GetLogz();
    TList primitives;
    TObjLink* link = pad->GetListOfPrimitives()->FirstLink();
    while (link) {
      TObject* o = link->GetObject();
      if (!o) { continue; }
      if (o->InheritsFrom("TH1") || o->InheritsFrom("TGraph")) {
        TObject* copy = o->Clone();
        Option_t* option = link->GetOption();
        primitives.Add(copy,option);
        name = o->GetName();
        title = o->GetTitle();
      }
      link = link->Next();
    }
    if (primitives.GetSize()==0) return;
    AliZMQMTviewerGUIview* window = new AliZMQMTviewerGUIview(
        name.c_str(),title.c_str(),100,200,700,600);
    window->fDrawnObjects.AddAll(&primitives);
    window->fCanvas.cd();
    window->fCanvas.SetLogx(logx);
    window->fCanvas.SetLogy(logy);
    window->fCanvas.SetLogz(logz);
    link = primitives.FirstLink();
    while (link) {
      TObject* o = link->GetObject();
      if (!o) { continue; }
      Option_t* option = link->GetOption();
      o->Draw(option);
      link = link->Next();
    }
    window->fCanvas.Update();
  }
}

//______________________________________________________________________________
void AliZMQMTviewerGUI::PadPicked(TPad* pad, TObject *object, Int_t event)
{
  if (pad) {
    fStatusBar->SetText(pad->GetName(),1);
  }
}

//______________________________________________________________________________
void AliZMQMTviewerGUI::ReconfigureViewer(std::string string)
{
  aliZMQmsg message;
  alizmq_msg_add(&message, "", string);
  alizmq_msg_send(&message, fZMQviewerConfig, 0);
  alizmq_msg_close(&message);
}

//______________________________________________________________________________
void AliZMQMTviewerGUI::DoPullButton()
{
  //just send something to interrupt sleep
  void* tmpsocket=NULL;
  int rc = 0;
  rc = alizmq_socket_init(tmpsocket, alizmq_context(), "PUSH>inproc://sleep", 0);
  rc = zmq_send(tmpsocket, &rc, sizeof(rc), 0);
  if (false) { printf("sleep interrupted rc: %i\n", rc); }
  zmq_close(tmpsocket);
}

//______________________________________________________________________________
void AliZMQMTviewerGUI::DoSelectionEntry()
{
  string sel = fSelectionEntry->GetText();
  string unsel = fUnSelectionEntry->GetText();
  delete fSelection; fSelection = NULL;
  delete fUnSelection; fUnSelection = NULL;
  if (sel.size()>0) fSelection = new TPRegexp(sel.c_str());
  if (unsel.size()>0) fUnSelection = new TPRegexp(unsel.c_str());

  if (fSelectionEntry->HasFocus()) {
    fUnSelectionEntry->SetFocus();
  } else if (fUnSelectionEntry->HasFocus()) {
    fSelectionEntry->SetFocus();
  }
  fViewer->GetSelection(&sel);
  fViewer->GetUnSelection(&unsel);
  bool tmp = true;
  fViewer->GetClearCanvas(&tmp);
  DoPullButton();
  UpdateCanvas();
}

ClassImp(AliZMQMTviewerGUIview)

void AliZMQMTviewerGUIview::CleanUp()
{
  fDrawnObjects.Delete();
  delete this;
}

