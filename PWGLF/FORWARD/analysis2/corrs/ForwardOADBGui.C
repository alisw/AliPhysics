
#ifndef __CINT__
#include "AliOADBForward.h"
#include <TGListBox.h>
#include <TGNumberEntry.h>
#include <TGTextEntry.h>
#include <TGComboBox.h>
#include <TGFrame.h>
#include <TGFileDialog.h>
#include <TGButtonGroup.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGMsgBox.h>
#include <TSystem.h>
#include <TError.h>
#include <TTimer.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TGListView.h>
#include <TGStatusBar.h>
#include <TDatime.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TBrowser.h>

namespace {
  void 
  ForwardOADBGUIErrorHandler(Int_t lvl, Bool_t doAbort, 
			     const char* location, 
			     const char* msg)
  {
    if (!doAbort && lvl >= kWarning) {
      EMsgBoxIcon msgIcon = kMBIconAsterisk;
      // if (lvl >= kInfo)    msgIcon = kMBIconAsterisk;
      if (lvl >= kWarning) msgIcon = kMBIconExclamation;
      if (lvl >= kError)   msgIcon = kMBIconStop;
      
      new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(), 
		   location, msg, msgIcon);
    }
    DefaultErrorHandler(lvl, doAbort, location, msg);
  }
}
#else 
class AliOADBForward;
class AliOADBForward::Entry;
// #if 0
class TGTransientFrame;
class TGFrame;
class TGLVEntry;
class TGHorizontalFrame;
class TGTextButton;
class TGCheckButton;
class TGTextEntry;
class TGVerticalFrame;
class TGLabel;
class TGComboBox;
class TGMainFrame;
class TGListView;
class TGLVContainer;
class TGHButtonGroup;
class TGLayoutHints;
class TGNumberEntry;
class TGStatusBar;
// #endif
#endif

struct ForwardOADBWait
{
  ForwardOADBWait(TGMainFrame* p)
    : fFrame(gClient->GetRoot(), p, 200, 40, kVerticalFrame),
      fHints(kLHintsExpandX,30,30,30,30),
      fLabel(&fFrame, ""), 
      // fProgress(&fFrame, 100), 
      fIsShown(false),
      fIncrement(0,false)
  {
    fFrame.AddFrame(&fLabel, &fHints );
    // fFrame.AddFrame(&fProgress, &fHints);
    fFrame.SetWindowName("Please wait ...");
    // fProgress.SetRange(0,1);
    fIncrement.Connect("Timeout()","ForwardOADBWait",this,"HandleIncr()");
  }
  ForwardOADBWait(const ForwardOADBWait&) {}
  ForwardOADBWait& operator=(const ForwardOADBWait&) { return *this; }
  void HandleIncr()
  {
    // Float_t dp = 0.1;
    // Float_t p  = fProgress.GetPosition();
    // Info("HandleIncr", "Handing increment (%f)", p);
    // if (p+dp >= 1) fProgress.SetPosition(0);
    // fProgress.Increment(dp);
    // fFrame.GetClient()->NeedRedraw(&fFrame);
  }
  void Show(const char* msg)
  {
    fLabel.SetText(msg);
    // fProgress.SetPosition(0);
    fFrame.MapSubwindows();
    Int_t width  = fFrame.GetDefaultWidth();
    Int_t height = fFrame.GetDefaultHeight();
    fFrame.Resize(width, height);
    fFrame.SetWMSize(width, height);
    fFrame.SetWMSizeHints(width, height, width, height, 0, 0);
    fFrame.MapRaised();
    fIsShown = true;
    // fIncrement.Start(0,false);
    fFrame.GetClient()->WaitForUnmap(&fFrame);    
  }
  void Hide() 
  {
    if (!fIsShown) return;
    // fIncrement.Stop();
    fFrame.UnmapWindow();
  }
  TGTransientFrame fFrame;
  TGLayoutHints    fHints;
  TGLabel          fLabel;
  Bool_t           fIsShown;
  TTimer           fIncrement;

};

struct ForwardOADBGUI
{
  enum { 
    kLabelWidth = 200
  };
  ForwardOADBGUI()
    : fMain(gClient->GetRoot(), 10, 10, kVerticalFrame),
      fOpenFrame(&fMain),
      fFileText(&fOpenFrame, "fmd_corrections.root"),
      fFileSelect(&fOpenFrame, "Browse"), 
      fTablesText(&fOpenFrame, "*"),
      fOpenButton(&fOpenFrame, "Open"),
      fRWButton(&fOpenFrame, "R/W"),
      fCloseButton(&fOpenFrame, "Close"),
      fSelectFrame(&fMain), 
      fTableFrame(&fSelectFrame), 
      fTableLabel(&fTableFrame, "Table: "), 
      fTableSelect(&fTableFrame),
      fRunFrame(&fSelectFrame), 
      fRunLabel(&fRunFrame, "Run: "), 
      fRunInput(&fRunFrame, 0, 0, -1, 
		TGNumberFormat::kNESReal,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMin, 0),
      fRunMode(&fRunFrame),
      fSysFrame(&fSelectFrame), 
      fSysLabel(&fSysFrame, "System: "), 
      fSysSelect(&fSysFrame),
      fSNNFrame(&fSelectFrame), 
      fSNNLabel(&fSNNFrame, "sqrt(sNN) [GeV]: "), 
      fSNNInput(&fSNNFrame, 0, 0, -1, TGNumberFormat::kNESReal,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMin, 0),
      fFldFrame(&fSelectFrame), 
      fFldLabel(&fFldFrame, "L3 field [kG]: "), 
      fFldSelect(&fFldFrame),
      fOtherFrame(&fSelectFrame),
      fMCButton(&fOtherFrame, "MC"),
      fSatButton(&fOtherFrame, "Satellite"),
      fOptionsFrame(&fSelectFrame), 
      fOptionsLabel(&fOptionsFrame, "Draw/Print options:"),
      fOptionsText(&fOptionsFrame, ""),
      fCommandFrame(&fSelectFrame), 
      fQueryButton(&fCommandFrame, "Query"),
      fListButton(&fCommandFrame, "List table"),
      fPrintButton(&fCommandFrame, "Print entry"),
      fCopyButton(&fCommandFrame, "Copy entry"),
      fDrawButton(&fCommandFrame, "Draw entry"),
      fPDFButton(&fCommandFrame, "Summarize entry"),
      fList(&fMain, 800, 400), 
      fListContainer(&fList),
      fFrameHints(kLHintsExpandX, 0, 0, 2, 0),
      fLabelHints(kLHintsNoHints, 4, 2, 0, 0),
      fEntryHints(kLHintsExpandX|kLHintsExpandY, 2, 4, 0, 0),
      fButtonHints(kLHintsExpandX, 2, 2, 0, 0),
      fListHints(kLHintsExpandX|kLHintsExpandY, 2, 2, 4, 2),
      fStatusBar(&fMain),
      fStatusBarHints(kLHintsExpandX, 2, 2, 4, 2),
      fMsg(&fMain),
      fDB(0),
      fEntry(0)
  {
    fMain.Connect("CloseWindow()", "ForwardOADBGUI", this, "HandleKill()");
    fMain.DontCallClose();

    fFileSelect.Connect("Clicked()", "ForwardOADBGUI", this, "HandleBrowse()");
    fOpenButton.Connect("Clicked()", "ForwardOADBGUI", this, "HandleOpen()");
    fCloseButton.Connect("Clicked()", "ForwardOADBGUI", this, "HandleClose()");
    fMain.AddFrame(&fOpenFrame, &fFrameHints);
    fOpenFrame.AddFrame(&fFileText, &fEntryHints);
    fOpenFrame.AddFrame(&fFileSelect, &fEntryHints);
    fOpenFrame.AddFrame(&fTablesText, &fEntryHints);
    fOpenFrame.AddFrame(&fOpenButton, &fEntryHints);
    fOpenFrame.AddFrame(&fRWButton,   &fButtonHints);
    fOpenFrame.AddFrame(&fCloseButton, &fEntryHints);

    fMain.AddFrame(&fSelectFrame, &fFrameHints);

    fTableLabel.SetWidth(kLabelWidth); fTableLabel.SetMinWidth(kLabelWidth);
    fTableSelect.SetHeight(22);
    fTableSelect.Connect("Selected(Int_t)","ForwardOADBGUI",this,
			 "HandleTable(Int_t)");
    fSelectFrame.AddFrame(&fTableFrame, &fFrameHints);
    fTableFrame.AddFrame(&fTableLabel, &fLabelHints);
    fTableFrame.AddFrame(&fTableSelect, &fEntryHints);
    
    fRunLabel.SetWidth(kLabelWidth); fRunLabel.SetMinWidth(kLabelWidth);
    fRunMode.AddEntry("default", 0);
    fRunMode.AddEntry("Exact",  1);
    fRunMode.AddEntry("Newest", 2);
    fRunMode.AddEntry("Near",   3);
    fRunMode.AddEntry("Older",  4);
    fRunMode.AddEntry("Newer",  5);
    fRunMode.SetHeight(22);
    fSelectFrame.AddFrame(&fRunFrame, &fFrameHints);
    fRunFrame.AddFrame(&fRunLabel, &fLabelHints);
    fRunFrame.AddFrame(&fRunInput, &fEntryHints);
    fRunFrame.AddFrame(&fRunMode, &fEntryHints);

    fSysLabel.SetWidth(kLabelWidth); fSysLabel.SetMinWidth(kLabelWidth);
    fSysSelect.AddEntry("- select -", 0);
    fSysSelect.AddEntry("p-p",   1);
    fSysSelect.AddEntry("Pb-Pb ",2);
    fSysSelect.AddEntry("p-Pb",  3);
    fSysSelect.AddEntry("Pb-p",  4);
    fSysSelect.AddEntry("Xe-Xe", 5);
    fSysSelect.SetHeight(22);
    fSelectFrame.AddFrame(&fSysFrame, &fFrameHints);
    fSysFrame.AddFrame(&fSysLabel, &fLabelHints);
    fSysFrame.AddFrame(&fSysSelect, &fEntryHints);

    fSNNLabel.SetWidth(kLabelWidth); fSNNLabel.SetMinWidth(kLabelWidth);
    fSNNInput.SetHeight(22);
    fSelectFrame.AddFrame(&fSNNFrame, &fFrameHints);
    fSNNFrame.AddFrame(&fSNNLabel, &fLabelHints);
    fSNNFrame.AddFrame(&fSNNInput, &fEntryHints);

    fFldLabel.SetWidth(kLabelWidth); fFldLabel.SetMinWidth(kLabelWidth);
    fFldSelect.AddEntry("- select -", 999);
    fFldSelect.AddEntry("-5", -5);
    fFldSelect.AddEntry("-2", -2);
    fFldSelect.AddEntry("0 ",  0);
    fFldSelect.AddEntry("+2", +2);
    fFldSelect.AddEntry("+5", +5);
    fFldSelect.SetHeight(22);
    fSelectFrame.AddFrame(&fFldFrame, &fFrameHints);
    fFldFrame.AddFrame(&fFldLabel, &fLabelHints);
    fFldFrame.AddFrame(&fFldSelect, &fEntryHints);

    fSelectFrame.AddFrame(&fOtherFrame, &fFrameHints);
    fOtherFrame.SetLayoutHints(&fButtonHints);
    fMCButton.AllowStayDown(true);
    fSatButton.AllowStayDown(true);
    // fOtherFrame.AddFrame(&fMCButton, &fEntryHints);
    // fOtherFrame.AddFrame(&fSatButton, &fEntryHints);
    // new TGCheckButton(&fOtherFrame, "MC:");
    // new TGCheckButton(&fOtherFrame, "Satellite:");

    fOptionsLabel.SetWidth(2*kLabelWidth);
    fSelectFrame.AddFrame(&fOptionsFrame, &fFrameHints);
    fOptionsFrame.AddFrame(&fOptionsLabel, &fLabelHints);
    fOptionsFrame.AddFrame(&fOptionsText, &fEntryHints);

    // Info("", "Connecting signals");
    fQueryButton.Connect("Clicked()", "ForwardOADBGUI", this, "HandleQuery()");
    fListButton.Connect("Clicked()", "ForwardOADBGUI", this, "HandleList()");
    fDrawButton.Connect("Clicked()", "ForwardOADBGUI", this, "HandleDraw()");
    fCopyButton.Connect("Clicked()", "ForwardOADBGUI", this, "HandleCopy()");
    fPDFButton.Connect("Clicked()", "ForwardOADBGUI", this, "HandlePDF()");
    fPrintButton.Connect("Clicked()", "ForwardOADBGUI", this, "HandlePrint()");
    fSelectFrame.AddFrame(&fCommandFrame, &fFrameHints);
    fCommandFrame.SetLayoutHints(&fButtonHints);
    
    // fList          = new TGListView(&fMain, 800, 400);
    // fListContainer = new TGLVContainer(fList);
    fListContainer.SetColHeaders("Entry", 
				 "Run", 
				 "System", 
				 "sqrt(sNN)", 
				 "L3 Field", 
				 "Type", 
				 "IP",
				 "Date",
				 "Author",
				 "AliROOT",
				 "Data");
    fList.SetViewMode(kLVDetails);
    fList.Connect("Clicked(TGLVEntry*,Int_t)", 
		   "ForwardOADBGUI", this, "HandleItem(TGLVEntry*,Int_t)");
    fList.Connect("DoubleClicked(TGLVEntry*,Int_t)", 
		   "ForwardOADBGUI", this, "HandleItem(TGLVEntry*,Int_t)");
    fListContainer.Connect("Clicked(TGFrame*,Int_t)",
			    "ForwardOADBGUI", this, 
			    "HandleItem(TGFrame*,Int_t)");
    fList.SetMinWidth(400);
    fList.SetMinHeight(200);
    fMain.AddFrame(&fList, &fListHints);

    fStatusBar.SetParts(1);
    fMain.AddFrame(&fStatusBar, &fStatusBarHints);
#ifndef __CINT__
    ::SetErrorHandler(ForwardOADBGUIErrorHandler);
#endif
    HandleEnable();

    fMain.SetMinWidth(600);
    fMain.SetMinHeight(480);
    fMain.MapSubwindows();    
    fMain.Resize(600, 480); // fMain.GetDefaultSize());
    fMain.MapWindow();
  }
  ~ForwardOADBGUI()
  {
    HandleClose();
#ifndef __CINT__
    ::SetErrorHandler(::DefaultErrorHandler);
#endif
    Info("~ForwardOADBGUI", "Closing");
  }
  ForwardOADBGUI(const ForwardOADBGUI&) 
    : fMain(gClient->GetRoot(), 10, 10, kVerticalFrame),
      fOpenFrame(&fMain),
      fFileText(&fOpenFrame, ""),
      fFileSelect(&fOpenFrame, ""), 
      fTablesText(&fOpenFrame, ""),
      fOpenButton(&fOpenFrame, ""),
      fRWButton(&fOpenFrame, ""),
      fCloseButton(&fOpenFrame, ""),
      fSelectFrame(&fMain), 
      fTableFrame(&fSelectFrame), 
      fTableLabel(&fTableFrame, ""), 
      fTableSelect(&fTableFrame),
      fRunFrame(&fSelectFrame), 
      fRunLabel(&fRunFrame, ""), 
      fRunInput(&fRunFrame),
      fRunMode(&fRunFrame),
      fSysFrame(&fSelectFrame), 
      fSysLabel(&fSysFrame, ""), 
      fSysSelect(&fSysFrame),
      fSNNFrame(&fSelectFrame), 
      fSNNLabel(&fSNNFrame, ""), 
      fSNNInput(&fSNNFrame),
      fFldFrame(&fSelectFrame), 
      fFldLabel(&fFldFrame, ""), 
      fFldSelect(&fFldFrame),
      fOtherFrame(&fSelectFrame),
      fMCButton(&fOtherFrame, ""),
      fSatButton(&fOtherFrame, ""),
      fOptionsFrame(&fSelectFrame), 
      fOptionsLabel(&fOptionsFrame, ""),
      fOptionsText(&fOptionsFrame, ""),
      fCommandFrame(&fSelectFrame), 
      fQueryButton(&fCommandFrame, ""),
      fListButton(&fCommandFrame, ""),
      fPrintButton(&fCommandFrame, ""),
      fCopyButton(&fCommandFrame, ""),
      fDrawButton(&fCommandFrame, ""),
      fPDFButton(&fCommandFrame, ""),
      fList(&fMain, 800, 400), 
      fListContainer(&fList),
      fFrameHints(),
      fLabelHints(),
      fEntryHints(),
      fButtonHints(),
      fListHints(),
      fStatusBar(&fMain),
      fStatusBarHints(),
      fMsg(0),
      fDB(0),
      fEntry(0)
  {}
  ForwardOADBGUI& operator=(const ForwardOADBGUI&) { return *this; }

  void UseDB(AliOADBForward* db)
  {
    if (!db) return;

    if (fDB) HandleClose();
    fEntry = 0;
    fDB = db;

    TString lt;
    const TMap& tables = fDB->GetTables();
    TIter next(&tables);
    TObject* key = 0;
    // Int_t    i   = 0;
    while ((key = next())) {
      AliOADBForward::Table* t = fDB->FindTable(key->GetName());
      
      lt.Append(Form("%s/%s", t->GetName(), t->fTree->GetTitle()));
    }
    fTablesText.SetText(lt);
    HandleEnable();    
  }
  void HandleKill()
  {
    // fMain.DontCallClose();
    fMain.DeleteWindow();
    Printf("Starting timer");
    // TTimer* t = new TTimer(Form("delete (ForwardOADBGUI*)%p", this), 100);
    // t->Start(100, true);
  }
  void HandleDBEntry(AliOADBForward::Entry* e)
  {
    Info("HandleDBEntry", "Selected entry %p", e);
    Bool_t en = (e != 0);
    fDrawButton.SetEnabled(en);
    fCopyButton.SetEnabled(en);
    fPrintButton.SetEnabled(en);
    fPDFButton.SetEnabled(en);
    
    fEntry = e;
  }
  void HandleEnable()
  {
    Bool_t enabled = fDB ? true : false;
    Bool_t hasTable = fTableSelect.GetSelected() != 0;

    fTableSelect.SetEnabled(enabled);
    fRunMode.SetEnabled(enabled);
    fSysSelect.SetEnabled(enabled);
    fFldSelect.SetEnabled(enabled);
    fMCButton.SetEnabled(enabled);
    fSatButton.SetEnabled(enabled);
    fQueryButton.SetEnabled(enabled && hasTable);
    fListButton.SetEnabled(enabled && hasTable);
    fPrintButton.SetEnabled(enabled && hasTable);
    fDrawButton.SetEnabled(enabled && hasTable);
    fCopyButton.SetEnabled(enabled && hasTable);
    fPDFButton.SetEnabled(enabled && hasTable);
    fOpenButton.SetEnabled(!enabled);
    fCloseButton.SetEnabled(enabled);
    HandleDBEntry(0);

    Int_t tsel = 0;
    if (!enabled) {
      fTableSelect.RemoveAll();
      fTableSelect.AddEntry("- select -", 0);
    }
    else {
      const TMap& tables = fDB->GetTables();
      TIter next(&tables);
      TObject* key = 0;
      Int_t    i   = 0;
      while ((key = next())) {
	fTableSelect.AddEntry(key->GetName(), ++i);
      }    
      if (tables.GetEntries() == 1) tsel = 1;
    }
    fTableSelect.Select(tsel, true);
    fSysSelect.Select(0, true);
    fFldSelect.Select(999, true);
    fRunMode.Select(0, true);

    fMain.Layout();
  }
  void HandleClose()
  {
    if (fDB) { 
      delete fDB;
      fDB = 0;
    }
    SetStatus("No DB connected");
    HandleEnable();
  }
  void HandleOpen()
  {
    if (fDB) HandleClose();
    Bool_t rw = fRWButton.IsOn();
    fDB = new AliOADBForward;
    Info("HandleOpen", "Opening DB file %s for tables %s", 
	 fFileText.GetText(), fTablesText.GetText());
    TString fn(fFileText.GetText());
    TString tn(fTablesText.GetText());
    if (!fDB->Open(fn, tn, rw, true, true)) { 
      Error("HandleOpen", "Failed to open database");
      delete fDB;
      fDB = 0;
    }
    SetStatus(Form("Connected to %s", fFileText.GetText()));
    // else 
    // fDB->Print();
    HandleEnable();
  }
  void HandleBrowse()
  {
    TGFileInfo fi;
    TString iniDir;
    if (gSystem->Getenv("OADB_PATH")) 
      iniDir = gSystem->ExpandPathName("$(OADB_PATH)");
    if (iniDir.IsNull()) 
      iniDir = gSystem->ExpandPathName("$(ALICE_PHYSICS)/OADB");
    iniDir.Append("/PWGLF/FORWARD/CORRECTIONS/data");
    char* ini = new char[iniDir.Length()+1];
    for (int i = 0; i < iniDir.Length(); i++) ini[i] = iniDir[i];
    ini[iniDir.Length()] = '\0';
    Printf("Initial directory: %s (%s)", iniDir.Data(), ini);
    fi.fIniDir = ini;
    new TGFileDialog(gClient->GetRoot(), &fMain, kFDOpen, &fi);

    TString nf = fi.fFilename; // 
    // nf = gSystem->ConcatFileName(fi.fIniDir, fi.fFilename);
    Info("HandleBrowse", "New file: %s", nf.Data());
    fFileText.SetText(nf);
  }
  void HandleEntry(Int_t i, AliOADBForward::Entry* e) 
  {
    TGLVEntry* lve = new TGLVEntry(&fListContainer, Form("%d", i), "");
    if (i < 0) lve->SetUserData(e);
    lve->SetUniqueID(i);
    TDatime dt(e->fTimestamp);
    lve->SetSubnames(Form("%lu", e->fRunNo), 
		     (e->fSys == 1 ? "p-p" : 
		      e->fSys == 2 ? "Pb-Pb" : 
		      e->fSys == 3 ? "p-Pb"  :
		      e->fSys == 4 ? "Pb-p"  :
		      e->fSys == 5 ? "Xe-Xe" :	"?"),
		     Form("%4huGeV",e->fSNN), 
		     Form("%+2hdkG", e->fField), 
		     (e->fMC ? "MC" : "Real"),
		     (e->fSatellite ? "Satellite" : "Nominal"), 
		     dt.AsSQLString(),
		     e->fAuthor, Form("%lu", e->fAliROOTRevision),
		     (e->fData ? e->fData->GetName() : "null"));
    fListContainer.AddItem(lve);
  }
  void HandleList()
  {
    if (!fDB) return;
    TString table;
    SelectedTable(table);

    if (table.IsNull()) {
      // Error("HandleList", "No table selected");
      return;
    }
    HandleDBEntry(0);
    AliOADBForward::Table* t= fDB->FindTable(table);
    if (!t) {
      Error("HandleList", "No table named %s in DB", table.Data());
      return;
    }
    // HandleQuery();
    t->Print(fOptionsText.GetText());
    // if (!fListContainer) return;
    
    fListContainer.RemoveAll();
    TTree* tree = t->fTree;
    Int_t  n    = tree->GetEntries();
    for (Int_t i = 0; i < n; i++) { 
      tree->GetEntry(i);
      AliOADBForward::Entry* e = t->fEntry;
      HandleEntry(i, e);
    }
    fList.AdjustHeaders();
    fMain.Layout();
  }
  void SelectedTable(TString& ret) const
  {
    ret = "";
    TGLBEntry* e = fTableSelect.GetSelectedEntry();
    if (!e) {
      Error("SelectedTable", "No table selected");
      return ;
    }
    ret = e->GetTitle();
  }
  void HandleTable(Int_t id)
  {
    Info("HandleTable", "Id=%d", id);
    fListButton.SetEnabled(id != 0);
    fQueryButton.SetEnabled(id != 0);
    fListContainer.RemoveAll();
    fList.AdjustHeaders();
    fMain.Layout();
    // fPrintButton.SetEnabled(enabled);
    // fDrawButton.SetEnabled(enabled);
    // fPDFButton.SetEnabled(enabled);
  }
  void HandleItem(TGFrame* lve, Int_t btn)
  {
    Info("HandleItem", "frame=%p", lve);
    HandleItem(static_cast<TGLVEntry*>(lve), btn);
  }
  void HandleItem(TGLVEntry* lve, Int_t)
  {
    Info("HandleItem", "entry=%p", lve);
    if (!lve) {
      Warning("HandleItem", "No item");
      return;
    }
    void* data = lve->GetUserData();
    AliOADBForward::Entry* e = 0;
    if (data) { 
      e = reinterpret_cast<AliOADBForward::Entry*>(data);
    }
    else { 
      TString tab;
      SelectedTable(tab);
      if (tab.IsNull()) return;
      
      AliOADBForward::Table* t = fDB->FindTable(tab);
      // Info("HandleItem", "Fetching item %d from table", lve->GetUniqueID());
      t->fTree->GetEntry(lve->GetUniqueID());
      e = t->fEntry;
    }
    if (!e) {
      Warning("HandleItem", "No entry");
      return;
    }
    // if (!gPad) TCanvas::MakeDefCanvas();
    e->Print();
    HandleDBEntry(e);
    // e->fData->Draw();
  }
  void HandlePrint() 
  {
    // TObject* o = HandleQuery();
    // if (!o) return;
    if (!fEntry) { 
      Warning("HandlePrint", "No entry selected");
      return;
    }
    fEntry->fData->Print(fOptionsText.GetText()); 
  }
  void HandleDraw()
  {
    if (!fEntry) { 
      Warning("HandleDraw", "No entry selected");
      return;
    }

    TObject* o = fEntry->fData;
    Info("HandleDraw", "Will draw object of type %s", o->ClassName());
    TString msg;
    if (o->IsA()->InheritsFrom(TParameter<Double_t>::Class())) { 
      TParameter<Double_t>* pd = static_cast<TParameter<Double_t>*>(o);
      msg = Form("%f", pd->GetVal());
    }
    if (o->IsA()->InheritsFrom(TParameter<Float_t>::Class())) { 
      TParameter<Float_t>* pf = static_cast<TParameter<Float_t>*>(o);
      msg = Form("%f", pf->GetVal());
    }
    if (!msg.IsNull()) {
      if (!gPad) TCanvas::MakeDefCanvas();
      
      TPaveText* t = new TPaveText(.1, .1, .9, .9);
      t->AddText(msg);
      t->Draw();
      return;
    }
    CorrDraw(o, false);
    // o->Draw(fOptionsText.GetText()); 
  }
  void HandleCopy()
  {
    TString table;
    SelectedTable(table);
    if (table.IsNull()) return;

    if (!fEntry) { 
      Warning("HandleCopy", "No entry selected");
      return;
    }

    TObject* o = fEntry->fData;
    Info("HandleCopy", "Will copy object of type %s", o->ClassName());


    ULong_t  newRun   = fRunInput.GetHexNumber();
    UShort_t newSys   = fSysSelect.GetSelected();
    UShort_t newSNN   = fSNNInput.GetIntNumber();
    Short_t  newFld   = fFldSelect.GetSelected();
    Bool_t   mc       = fMCButton.IsDown();
    Bool_t   sat      = fSatButton.IsDown();
    ULong_t  oldRun   = fEntry->fRunNo;
    UShort_t oldSys   = fEntry->fSys;
    UShort_t oldSNN   = fEntry->fSNN;
    Short_t  oldFld   = fEntry->fField;

    TString msg;
    msg = Form("Will copy %lu/%hu/%hu/%hd -> %lu/%hu/%hu/%hd (%s,%s) in %s",
	       oldRun, oldSys, oldSNN, oldFld, 
	       newRun, newSys, newSNN, newFld, 
	       (mc ? "mc" : "real"), (sat ? "satellite" : "nominal"),
	       table.Data());
    
    Int_t ret = 0;
    new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(), 
		 "Confirm copy of OADB entry",
		 msg.Data(), 0, kMBOk|kMBCancel, &ret);
    
    if (ret ==  kMBCancel) {
      Info("HandleCopy", "%s CANCELLED", msg.Data());
      return;
    }
    Info("HandleCopy", "%s CONFIRMED", msg.Data());

    if (!fDB->CopyEntry(table, 
			oldRun, oldSys, oldSNN, oldFld, 
			newRun, newSys, newSNN, newFld, 
			mc, sat)) { 
      // Warning("HandleCopy", "Copy failed!");
      return;
    }
    fDB->Update();
    

    // o->Draw(fOptionsText.GetText()); 
  }
  void HandlePDF()
  {
    if (!fEntry) { 
      Warning("HandlePrint", "No entry selected");
      return;
    }


    TObject* o = fEntry->fData;
    gROOT->SetBatch(true);
    CorrDraw(o, true);
    gROOT->SetBatch(false);
    // fEntry->fData->SaveAs(out, fOptionsText.GetText()); 

  }
  void LoadCorrDraw()
  {
    const char* opt = "++g";
    const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
    gSystem->AddIncludePath(Form("-I$ALICE_ROOT/include "
				 "-I$ALICE_PHYSICS/include "
				 "-I%s -I%s/scripts", fwd, fwd));
    Info("CorrDraw", "Loading SummaryDrawer.C%s", opt);
    gROOT->LoadMacro(Form("%s/scripts/SummaryDrawer.C%s", fwd, opt));
    // gROOT->ProcessLine(".Class SummaryDrawer");
    Info("CorrDraw", "Loading CorrDrawer.C%s", opt);
    gROOT->LoadMacro(Form("%s/corrs/CorrDrawer.C%s", fwd, opt));
    // gROOT->ProcessLine(".Class SummaryDrawer");
    // gROOT->ProcessLine(".Class CorrDrawer");
    CloseMsg();
  }    
  void CorrDraw(const TObject* o, Bool_t summarize)
  {
    if (!gROOT->GetClass("CorrDrawer")) { 
      TTimer* t = new TTimer(Form("((ForwardOADBGUI*)%p)->LoadCorrDraw()",this),
			     0, kTRUE);
      t->Start(100, true);
      MakeMsg("Compiling drawer");
    }
    o->Print();
    if (summarize) {
      TTimer* t = new TTimer(Form("((ForwardOADBGUI*)%p)"
				  "->DoCorrDraw((TObject*)%p,1)",this, o), 
			     0, kTRUE);
      t->Start(100, true);
      MakeMsg("Drawing to PDF");
    }
    else 
      DoCorrDraw(o, false);
  }
  void DoCorrDraw(const TObject* o, Bool_t summarize) 
  {
    TString cmd(Form("CorrDrawer cd; cd.Summarize((const %s*)%p,%d);",
		     o->ClassName(), o, summarize));
    gROOT->ProcessLine(cmd);
    CloseMsg();
  }  
  void MakeFileName(TString& out) const
  {
    if (!fEntry) return;
    
    SelectedTable(out);
    if (out.IsNull()) return;
    
    out.Append(Form("_run%09lu", fEntry->fRunNo));
    out.Append(Form("_%s", (fEntry->fSys == 1 ? "pp" : 
			    fEntry->fSys == 2 ? "PbPb" :
			    fEntry->fSys == 3 ? "pPb" :
			    fEntry->fSys == 4 ? "Pbp" :
			    fEntry->fSys == 5 ? "Xe-Xe" : "XX")));
    out.Append(Form("_%04huGeV", fEntry->fSNN));
    out.Append(Form("_%c%hukG", fEntry->fField >= 0 ? 'p' : 'm', 
		    TMath::Abs(fEntry->fField)));
    out.Append(Form("_%s", fEntry->fMC ? "mc" : "real"));
    out.Append(Form("_%s", fEntry->fSatellite ? "sat" : "nom"));
    out.Append(".pdf");
  }
  TObject* HandleQuery()
  {
    ULong_t  run   = fRunInput.GetHexNumber();
    Short_t  mode  = fRunMode.GetSelected();
    Short_t  sys   = fSysSelect.GetSelected();
    UShort_t sNN   = fSNNInput.GetIntNumber();
    Short_t  fld   = fFldSelect.GetSelected();
    Bool_t   mc    = fMCButton.IsDown();
    Bool_t   sat   = fSatButton.IsDown();
    TString  tab;
    SelectedTable(tab);
    
    Info("HandleQuery", "tab=%s runNo=%lu mode=%d sys=%d "
	 "sNN=%d fld=%d mc=%d sat=%d", 
	 tab.Data(), run, mode, sys, sNN, fld, mc, sat);

    if (tab.IsNull()) {
      // Error("HandleQuery", "No table selected");
      return 0;
    }
    AliOADBForward::ERunSelectMode qmode = AliOADBForward::kDefault;
    switch (mode) { 
    case 0: qmode = AliOADBForward::kDefault; break;
    case 1: qmode = AliOADBForward::kExact; break;
    case 2: qmode = AliOADBForward::kNewest; break;
    case 3: qmode = AliOADBForward::kNear; break;
    case 4: qmode = AliOADBForward::kOlder; break;
    case 5: qmode = AliOADBForward::kNewer; break;
    }

    Info("HandleQuery", "tab=%s runNo=%lu mode=%d sys=%d "
	 "sNN=%d fld=%d mc=%d sat=%d", 
	 tab.Data(), run, qmode, sys, sNN, fld, mc, sat);

    AliOADBForward::Entry* e = fDB->Get(tab, run, qmode, sys, sNN, 
					fld, mc, sat);
    if (!e) return 0;
    // if (drawNotPrint) e->Inspect();
    // else              e->Print(fOptionsText.GetText());
    e->Print();
    // if (fListContainer) { 
    fListContainer.RemoveAll();
    HandleEntry(-1, e);
    // }
    if (!e->fData) return 0;
    HandleDBEntry(e);

    fList.AdjustHeaders();
    fMain.Layout();

    return e->fData;
  }
  TGMainFrame* GetMain() { return &fMain; }
  void MakeMsg(const char* what) 
  { 
    fMsg.Show(what);
  }
  void CloseMsg()
  {
    Info("CloseMsg", "Closing message window");
    fMsg.Hide();
  }
  void SetStatus(const char* what)
  {
    fStatusBar.SetText(what,0);
  }
  TGMainFrame       fMain;
  TGHorizontalFrame fOpenFrame;
  TGTextEntry       fFileText;
  TGTextButton      fFileSelect;
  TGTextEntry       fTablesText;
  TGTextButton      fOpenButton;
  TGCheckButton     fRWButton;
  TGTextButton      fCloseButton;
  TGVerticalFrame   fSelectFrame;
  TGHorizontalFrame fTableFrame;
  TGLabel           fTableLabel;
  TGComboBox        fTableSelect;
  TGHorizontalFrame fRunFrame;
  TGLabel           fRunLabel;
  TGNumberEntry     fRunInput;
  TGComboBox        fRunMode;
  TGHorizontalFrame fSysFrame;
  TGLabel           fSysLabel;
  TGComboBox        fSysSelect;
  TGHorizontalFrame fSNNFrame;
  TGLabel           fSNNLabel;
  TGNumberEntry     fSNNInput;
  TGHorizontalFrame fFldFrame;
  TGLabel           fFldLabel;
  TGComboBox        fFldSelect;
  TGHButtonGroup    fOtherFrame;
  TGTextButton      fMCButton;
  TGTextButton      fSatButton;
  TGHorizontalFrame fOptionsFrame;
  TGLabel           fOptionsLabel;
  TGTextEntry       fOptionsText;
  TGHButtonGroup    fCommandFrame;
  TGTextButton      fQueryButton;  
  TGTextButton      fListButton;
  TGTextButton      fPrintButton;
  TGTextButton      fCopyButton;
  TGTextButton      fDrawButton;
  TGTextButton      fPDFButton;
  TGListView        fList;
  TGLVContainer     fListContainer;
  TGLayoutHints     fFrameHints;
  TGLayoutHints     fLabelHints;
  TGLayoutHints     fEntryHints; 
  TGLayoutHints     fButtonHints;
  TGLayoutHints     fListHints;
  TGStatusBar       fStatusBar;
  TGLayoutHints     fStatusBarHints;
  ForwardOADBWait   fMsg;
  AliOADBForward*   fDB;
  AliOADBForward::Entry* fEntry;
  // TCanvas*          fDataCanvas;
  // TCanvas*          fEntryCanvas;
};


TGMainFrame* ForwardOADBGui(AliOADBForward* db=0,
			    const char* path=0)
{
  const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
  if (!gROOT->GetClass("AliOADBForward")) 
    // gSystem->Load("libGui");
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
  
  // gSystem->AddIncludePath(Form("-I%s", fwd));
  // gROOT->LoadMacro(Form("%s/corrs/ForwardOADBGUI.C", fwd));

  new TBrowser;
  // new TGClient();
  ForwardOADBGUI* gui = new ForwardOADBGUI();
  if (db) {
    gui->UseDB(db);
    if (path) gui->SetStatus(Form("Connected to %s", path));
  }
  return gui->GetMain();
}
  
//
// EOF
//
  
