
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TBrowser.h"
#include "TGClient.h"
#include "TGFrame.h"
#include "TGTab.h"
#include "TGMenu.h"
#include "TGLayout.h"
#include "TGSplitter.h"
#include "TGStatusBar.h"
#include "Varargs.h"
#include "TInterpreter.h"
#include "TBrowser.h"
#include "TGFileDialog.h"
#include "TObjString.h"
#include <KeySymbols.h>

#include "TGNewBrowser.h"
#include "TGFileBrowser.h"

#include "Getline.h"

#ifdef WIN32
#include <TWin32SplashThread.h>
#endif

static const char *gOpenFileTypes[] = {
   "ROOT files",   "*.root",
   "All files",    "*",
   0,              0 
};

static const char *gPluginFileTypes[] = {
   "ROOT files",   "*.C",
   "All files",    "*",
   0,              0 
};

enum ENewBrowserMessages {
   kBrowse = 11011,
   kOpenFile,
   kNewEditor,
   kNewCanvas,
   kExecPlugin,
   kCloseTab,
   kCloseWindow,
   kQuitRoot
};

//______________________________________________________________________________
TGNewBrowser::TGNewBrowser(const char *name, UInt_t width, UInt_t height, Bool_t initshow)
   : TGMainFrame(gClient->GetDefaultRoot(), width, height)
{
   // Create browser with a specified width and height.

   CreateBrowser(name);
   Resize(width, height);
   if (initshow) {
      InitPlugins();
      MapWindow();
   }
   gVirtualX->SetInputFocus(GetId());
}

//______________________________________________________________________________
TGNewBrowser::TGNewBrowser(const char *name, Int_t x, Int_t y,
                           UInt_t width, UInt_t height, Bool_t initshow)
   : TGMainFrame(gClient->GetDefaultRoot(), width, height)
{
   // Create browser with a specified width and height and at position x, y.

   CreateBrowser(name);
   MoveResize(x, y, width, height);
   SetWMPosition(x, y);
   if (initshow) {
      InitPlugins();
      MapWindow();
   }
   gVirtualX->SetInputFocus(GetId());
}


//______________________________________________________________________________
void TGNewBrowser::CreateBrowser(const char *name)
{

   fVf = new TGVerticalFrame(this, 100, 100);

   fLH0 = new TGLayoutHints(kLHintsNormal);
   fLH1 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fLH2 = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 1, 1, 1, 3);
   fLH3 = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX);
   fLH4 = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2);
   fLH5 = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY);
   fLH6 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX);
   fLH7 = new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY);

   // Menubar Frame
   fTopMenuFrame = new TGHorizontalFrame(fVf, 100, 20);

   fPreMenuFrame = new TGHorizontalFrame(fTopMenuFrame, 0, 20, kRaisedFrame);
   fMenuBar   = new TGMenuBar(fPreMenuFrame, 10, 10, kHorizontalFrame);
   fMenuFile  = new TGPopupMenu(gClient->GetDefaultRoot());
   fMenuFile->AddEntry("&Browse...             Ctrl+B", kBrowse);
   fMenuFile->AddEntry("&Open...                 Ctrl+O", kOpenFile);
   fMenuFile->AddEntry("New &Editor          Ctrl+E", kNewEditor);
   fMenuFile->AddEntry("New &Canvas       Ctrl+C", kNewCanvas);
   fMenuFile->AddEntry("Exec &Plugin...      Ctrl+P", kExecPlugin);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("Close &Tab           Ctrl+T", kCloseTab);
   fMenuFile->AddEntry("Close &Window   Ctrl+W", kCloseWindow);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("&Quit Root             Ctrl+Q", kQuitRoot);
   fMenuBar->AddPopup("Framewor&k", fMenuFile, fLH1);
   fMenuFile->Connect("Activated(Int_t)", "TGNewBrowser", this,
                      "HandleMenu(Int_t)");
   fPreMenuFrame->AddFrame(fMenuBar, fLH2);
   fTopMenuFrame->AddFrame(fPreMenuFrame, fLH0);

   fMenuFrame = new TGHorizontalFrame(fTopMenuFrame, 100, 20, kRaisedFrame);
   fTopMenuFrame->AddFrame(fMenuFrame, fLH5);

   fVf->AddFrame(fTopMenuFrame, fLH3);
   fActMenuBar = fMenuBar;

   // Toolbar Frame
   fToolbarFrame = new TGHorizontalFrame(fVf, 100, 20, kHorizontalFrame | 
                                         kRaisedFrame);
   fVf->AddFrame(fToolbarFrame, fLH3);

   fHf = new TGHorizontalFrame(fVf, 100, 100);
   // Tabs & co...
   fV1 = new TGVerticalFrame(fHf, 250, 100, kFixedWidth);
   fV2 = new TGVerticalFrame(fHf, 600, 100);
   fH1 = new TGHorizontalFrame(fV2, 100, 100);
   fH2 = new TGHorizontalFrame(fV2, 100, 100, kFixedHeight);

   // Left tab
   fTabLeft = new TGTab(fV1,100,100);
   //fTabLeft->AddTab("Tab 1");
   fTabLeft->Resize(fTabLeft->GetDefaultSize());
   fV1->AddFrame(fTabLeft, fLH4);

   // Vertical splitter   
   fVSplitter = new TGVSplitter(fHf, 4, 4);
   fVSplitter->SetFrame(fV1, kTRUE);
   fHf->AddFrame(fV1, fLH7);
   fHf->AddFrame(fVSplitter, fLH7);

   // Right tab
   fTabRight = new TGTab(fH1, 500, 100);
   //fTabRight->AddTab("Tab 1");
   fTabRight->Resize(fTabRight->GetDefaultSize());
   fH1->AddFrame(fTabRight, fLH5);
   fTabRight->Connect("Selected(Int_t)", "TGNewBrowser", this, "DoTab(Int_t)");
   fV2->AddFrame(fH1, fLH4);
   
   // Horizontal splitter
   fHSplitter = new TGHSplitter(fV2, 4, 4);
   fV2->AddFrame(fHSplitter, fLH3);

   // Bottom tab
   fTabBottom = new TGTab(fH2, 100, 100);
   //fTabBottom->AddTab("Tab 1");
   fH2->AddFrame(fTabBottom, fLH4);
   fV2->AddFrame(fH2, fLH3);
   
   fHSplitter->SetFrame(fH2, kFALSE);
   fHf->AddFrame(fV2, fLH5);
   fVf->AddFrame(fHf, fLH5);
   AddFrame(fVf, fLH5);

   // status bar
   fStatusBar = new TGStatusBar(this, 400, 20);
   Int_t parts[] = { 80, 20 };
   fStatusBar->SetParts(parts, 2);
   AddFrame(fStatusBar, fLH6);

   fEditFrame = 0;
   fEditTab   = 0;
   fNbTab[0]  = fNbTab[1] = fNbTab[2] = 0;
   fCrTab[0]  = fCrTab[1] = fCrTab[2] = -1;
   
   // Set a name to the main frame   
   SetWindowName(name);
   SetIconName(name);
   SetClassHints("Browser", "Browser");

   SetWMSizeHints(600, 350, 10000, 10000, 2, 2);
   MapSubwindows();
   Resize(GetDefaultSize());
   AddInput(kKeyPressMask | kKeyReleaseMask);

   fVf->HideFrame(fToolbarFrame);
}  

//______________________________________________________________________________
TGNewBrowser::~TGNewBrowser()
{
   // Clean up all widgets, frames and layouthints that were used

   delete fLH0;
   delete fLH1;
   delete fLH2;
   delete fLH3;
   delete fLH4;
   delete fLH5;
   delete fLH6;
   delete fLH7;
   delete fMenuFile;
   delete fMenuBar;
   delete fMenuFrame;
   delete fToolbarFrame;
   delete fVSplitter;
   delete fHSplitter;
   delete fTabLeft;
   delete fTabRight;
   delete fTabBottom;
   delete fH1;
   delete fH2;
   delete fV1;
   delete fV2;
   delete fHf;
   delete fStatusBar;
   delete fVf;
}

//______________________________________________________________________________
void TGNewBrowser::CloseWindow()
{
   // Called when window is closed via the window manager.

   TGFrameElement *el;
   Int_t i;
   Disconnect(fMenuFile, "Activated(Int_t)", this, "HandleMenu(Int_t)");
   for (i=0;i<fTabLeft->GetNumberOfTabs();i++) {
      el = (TGFrameElement *)fTabLeft->GetTabContainer(i)->GetList()->First();
      if (el && el->fFrame) {
         el->fFrame->SetFrameElement(0);
         if (el->fFrame->InheritsFrom("TGMainFrame"))
            ((TGMainFrame *)el->fFrame)->CloseWindow();
         else
            delete el->fFrame;
         el->fFrame = 0;
         if (el->fLayout && (el->fLayout != fgDefaultHints) &&
            (el->fLayout->References() > 0)) {
            el->fLayout->RemoveReference();
         }
         fTabLeft->GetTabContainer(i)->GetList()->Remove(el);
         delete el;
      }
   }
   for (i=0;i<fTabRight->GetNumberOfTabs();i++) {
      el = (TGFrameElement *)fTabRight->GetTabContainer(i)->GetList()->First();
      if (el && el->fFrame) {
         el->fFrame->SetFrameElement(0);
         if (el->fFrame->InheritsFrom("TGMainFrame"))
            ((TGMainFrame *)el->fFrame)->CloseWindow();
         else
            delete el->fFrame;
         el->fFrame = 0;
         if (el->fLayout && (el->fLayout != fgDefaultHints) &&
            (el->fLayout->References() > 0)) {
            el->fLayout->RemoveReference();
         }
         fTabRight->GetTabContainer(i)->GetList()->Remove(el);
         delete el;
      }
   }
   for (i=0;i<fTabBottom->GetNumberOfTabs();i++) {
      el = (TGFrameElement *)fTabBottom->GetTabContainer(i)->GetList()->First();
      if (el && el->fFrame) {
         el->fFrame->SetFrameElement(0);
         if (el->fFrame->InheritsFrom("TGMainFrame"))
            ((TGMainFrame *)el->fFrame)->CloseWindow();
         else
            delete el->fFrame;
         el->fFrame = 0;
         if (el->fLayout && (el->fLayout != fgDefaultHints) &&
            (el->fLayout->References() > 0)) {
            el->fLayout->RemoveReference();
         }
         fTabBottom->GetTabContainer(i)->GetList()->Remove(el);
         delete el;
      }
   }
   DeleteWindow();
}

//______________________________________________________________________________
void TGNewBrowser::DoTab(Int_t id)
{
   // Handle Tab navigation.

   TGTab *sender = (TGTab *)gTQSender;
   if ((sender) && (sender == fTabRight)) {
      SwitchMenus(sender->GetTabContainer(id));
   }
}

//______________________________________________________________________________
void TGNewBrowser::ExecPlugin(const char *fname, Int_t pos, Int_t subpos)
{
   // Execute a macro and embed the created frame in the tab "pos" 
   // and tab element "subpos".

   StartEmbedding(pos, subpos);
   gROOT->Macro(fname);
   StopEmbedding();
}

//______________________________________________________________________________
Bool_t TGNewBrowser::HandleKey(Event_t *event)
{
   char   input[10];
   Int_t  n;
   UInt_t keysym;

   if (event->fType == kGKeyPress) {
      gVirtualX->LookupString(event, input, sizeof(input), keysym);
      n = strlen(input);

      switch ((EKeySym)keysym) {   // ignore these keys
         case kKey_Shift:
         case kKey_Control:
         case kKey_Meta:
         case kKey_Alt:
         case kKey_CapsLock:
         case kKey_NumLock:
         case kKey_ScrollLock:
            return kTRUE;
         default:
            break;
      }
      if (event->fState & kKeyControlMask) {   // Cntrl key modifier pressed
         switch ((EKeySym)keysym & ~0x20) {   // treat upper and lower the same
            case kKey_B:
               fMenuFile->Activated(kBrowse);
               return kTRUE;
            case kKey_O:
               fMenuFile->Activated(kOpenFile);
               return kTRUE;
            case kKey_E:
               fMenuFile->Activated(kNewEditor);
               return kTRUE;
            case kKey_C:
               fMenuFile->Activated(kNewCanvas);
               return kTRUE;
            case kKey_P:
               fMenuFile->Activated(kExecPlugin);
               return kTRUE;
            case kKey_T:
               fMenuFile->Activated(kCloseTab);
               return kTRUE;
            case kKey_W:
               fMenuFile->Activated(kCloseWindow);
               return kTRUE;
            case kKey_Q:
               fMenuFile->Activated(kQuitRoot);
               return kTRUE;
            default:
               break;
         }
      }
   }
   return TGMainFrame::HandleKey(event);
}

//______________________________________________________________________________
void TGNewBrowser::HandleMenu(Int_t id)
{
   // Handle menu items.

   static Int_t cNr = 1;
   static Int_t eNr = 1;
   TGPopupMenu *sender = (TGPopupMenu *)gTQSender;
   if (sender != fMenuFile)
      return;
   switch (id) {
      case kBrowse:
         new TBrowser();
         break;
      case kOpenFile:
         {
            static TString dir(".");
            TGFileInfo fi;
            fi.fFileTypes = gOpenFileTypes;
            fi.fIniDir    = StrDup(dir);
            new TGFileDialog(gClient->GetDefaultRoot(), this,
                             kFDOpen,&fi);
            dir = fi.fIniDir;
            if (fi.fMultipleSelection && fi.fFileNamesList) {
               TObjString *el;
               TIter next(fi.fFileNamesList);
               while ((el = (TObjString *) next())) {
                  gROOT->ProcessLine(Form("new TFile(\"%s\");",
                                     gSystem->UnixPathName(el->GetString())));
               }
            }
            else if (fi.fFilename) {
               gROOT->ProcessLine(Form("new TFile(\"%s\");",
                                  gSystem->UnixPathName(fi.fFilename)));
            }
         }
         break;
      case kNewEditor:
         StartEmbedding(1);
         ++eNr;
         gROOT->ProcessLine(Form("new TGTextEditor((const char *)0, (const TGWindow *)0x%lx)", 
                            gClient->GetRoot()));
         StopEmbedding();
         SetTabTitle(Form("Editor %d", eNr), 1);
         break;
      case kNewCanvas:
         StartEmbedding(1);
         ++cNr;
         gROOT->ProcessLine(Form("new TCanvas(\"BrowserCanvas%d\", \"Browser Canvas %d\")", cNr, cNr));
         StopEmbedding();
         SetTabTitle(Form("Canvas %d", cNr), 1);
         break;
      case kExecPlugin:
         {
            static TString dir(".");
            TGFileInfo fi;
            fi.fFileTypes = gPluginFileTypes;
            fi.fIniDir    = StrDup(dir);
            new TGFileDialog(gClient->GetDefaultRoot(), this,
                             kFDOpen,&fi);
            dir = fi.fIniDir;
            if (fi.fFilename) {
               ExecPlugin(fi.fFilename, kRight);
            }
         }
         break;
      case kCloseTab:
         RemoveFrame(kRight, fTabRight->GetCurrent());
         break;
      case kCloseWindow:
         CloseWindow();
         break;
      case kQuitRoot:
         gApplication->Terminate(0);
         break;
      default:
         break;
   }
}

//______________________________________________________________________________
void TGNewBrowser::RemoveFrame(Int_t pos, Int_t subpos)
{
   // Remove tab element "subpos" from tab "pos".

   TGTab *edit = 0;
   switch (pos) {
      case kLeft: // left
         edit = fTabLeft;
         break;
      case kRight: // right
         edit = fTabRight;
         fMenuFrame->HideFrame(fActMenuBar);
         fMenuFrame->GetList()->Remove(fActMenuBar);
         fActMenuBar = 0;
         break;
      case kBottom: // bottom
         edit = fTabBottom;
         break;
   }
   TGFrameElement *el = (TGFrameElement *)edit->GetTabContainer(subpos)->GetList()->First();
   if (el && el->fFrame) {
      el->fFrame->SetFrameElement(0);
      if (el->fFrame->InheritsFrom("TGMainFrame"))
         ((TGMainFrame *)el->fFrame)->CloseWindow();
      else
         delete el->fFrame;
      el->fFrame = 0;
      if (el->fLayout && (el->fLayout != fgDefaultHints) &&
         (el->fLayout->References() > 0)) {
         el->fLayout->RemoveReference();
      }
      edit->GetTabContainer(subpos)->GetList()->Remove(el);
      delete el;
   }
   fNbTab[pos]--;
   edit->RemoveTab(subpos);
}

//______________________________________________________________________________
TGTab* TGNewBrowser::GetTab(Int_t pos) const
{
   switch (pos) {
      case kLeft:   return fTabLeft;
      case kRight:  return fTabRight;
      case kBottom: return fTabBottom;
      default:      return 0;
   }
}

//______________________________________________________________________________
void TGNewBrowser::SetTab(Int_t pos, Int_t subpos)
{
   TGTab *tab = GetTab(pos);
   if (subpos == -1)
      subpos = fCrTab[pos];

   if (tab->SetTab(subpos, kFALSE)) { // Block signal emit
      if (pos == kRight)
         SwitchMenus(tab->GetTabContainer(subpos));
      tab->Layout();
   }
}

//______________________________________________________________________________
void TGNewBrowser::SetTabTitle(const char *title, Int_t pos, Int_t subpos)
{
   TGTab *edit = GetTab(pos);
   if (subpos == -1)
      subpos = fCrTab[pos];

   TGTabElement *el = edit->GetTabTab(subpos);
   if (el) {
      el->SetText(new TGString(title));
      edit->Layout();
   }
}

//______________________________________________________________________________
void TGNewBrowser::ShowMenu(TGCompositeFrame *menu)
{
   // Show the selected frame's menu and hide previous one.

   TGFrameElement *el = 0;
   // temporary solution until I find a proper way to handle 
   // these bloody menus...
   fBindList->Delete();
   TIter nextm(fMenuBar->GetList());
   while ((el = (TGFrameElement *) nextm())) {
      TGMenuTitle *t = (TGMenuTitle *) el->fFrame;
      Int_t code = t->GetHotKeyCode();
      BindKey(fMenuBar, code, kKeyMod1Mask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyShiftMask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyLockMask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyShiftMask | kKeyLockMask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyMod2Mask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyMod2Mask | kKeyLockMask);
      BindKey(fMenuBar, code, kKeyMod1Mask | kKeyShiftMask | kKeyMod2Mask | kKeyLockMask);
   }
   fMenuFrame->HideFrame(fActMenuBar);
   fMenuFrame->ShowFrame(menu);
   menu->Layout();
   fMenuFrame->Layout();
   fActMenuBar = menu;
}

//______________________________________________________________________________
TGCompositeFrame* TGNewBrowser::StartEmbedding(Int_t pos, Int_t subpos)
{
   // Start embedding external frame in the tab "pos" and tab element "subpos".

   fEditTab = GetTab(pos);

   if (fEditFrame == 0) {
      if (subpos == -1) {
         fCrTab[pos] = fNbTab[pos]++;
         fEditFrame = fEditTab->AddTab(Form("Tab %d",fNbTab[pos]));
         fEditFrame->MapWindow();
         TGTabElement *tabel = fEditTab->GetTabTab(fEditTab->GetNumberOfTabs()-1);
         if(tabel) tabel->MapWindow();
         fEditTab->SetTab(fEditTab->GetNumberOfTabs()-1);
         fEditTab->Layout();
      }
      else {
         fCrTab[pos] = subpos;
         fEditFrame = fEditTab->GetTabContainer(subpos);
         fEditTab->SetTab(subpos);
      }
      fEditFrame->SetEditable();
   }
   return fEditFrame;
}

//______________________________________________________________________________
void TGNewBrowser::StopEmbedding(TGLayoutHints *layout)
{
   // Stop embedding external frame in the current editable frame.
   
   if (fEditFrame != 0) {
      fEditFrame->SetEditable(kFALSE);
      if (layout) {
         TGFrameElement *el = (TGFrameElement*) fEditFrame->GetList()->Last();
         // !!!! MT what to do with the old layout? Leak it for now ...
         el->fLayout = layout;
      }
      fEditFrame->Layout();
      if (fEditTab == fTabRight)
         SwitchMenus(fEditFrame);
      fEditFrame = fEditTab = 0;
   }
}

//______________________________________________________________________________
void TGNewBrowser::SwitchMenus(TGCompositeFrame  *from)
{
   // Move the menu from original frame to our TGMenuFrame, or display the 
   // menu associated to the current tab.

   TGFrameElement *fe = (TGFrameElement *)from->GetList()->First();
   if (!fe) {
      if (fActMenuBar != fMenuBar)
         ShowMenu(fMenuBar);
      return;
   }
   TGCompositeFrame *embed = (TGCompositeFrame *)fe->fFrame;
   TGFrameElement *el = 0;
   if (embed && embed->GetList()) {
      TIter next(embed->GetList());
      while ((el = (TGFrameElement *)next())) {
         if (el->fFrame->InheritsFrom("TGMenuBar")) {
            TGMenuBar *menu = (TGMenuBar *)el->fFrame;
            if (fActMenuBar == menu)
               return;
            TGFrameElement *nw;
            TIter nel(fMenuFrame->GetList());
            while ((nw = (TGFrameElement *) nel())) {
               if (nw->fFrame == menu) {
                  ShowMenu(menu);
                  return;
               }
            }
            ((TGCompositeFrame *)menu->GetParent())->HideFrame(menu);
            menu->ReparentWindow(fMenuFrame);
            fMenuFrame->AddFrame(menu, fLH2);
            TGFrameElement *mel;
            TIter mnext(menu->GetList());
            while ((mel = (TGFrameElement *) mnext())) {
               TGMenuTitle *t = (TGMenuTitle *) mel->fFrame;
               TGPopupMenu *popup = menu->GetPopup(t->GetName());
               RecursiveReparent(popup);
               if (popup->GetEntry("Close Canvas")) {
                  TGMenuEntry *exit = popup->GetEntry("Close Canvas");
                  popup->HideEntry(exit->GetEntryId());
               }
               if (popup->GetEntry("Close Viewer")) {
                  TGMenuEntry *exit = popup->GetEntry("Close Viewer");
                  popup->HideEntry(exit->GetEntryId());
               }
               if (popup->GetEntry("Quit ROOT")) {
                  TGMenuEntry *exit = popup->GetEntry("Quit ROOT");
                  popup->HideEntry(exit->GetEntryId());
               }
               if (popup->GetEntry("Exit")) {
                  TGMenuEntry *exit = popup->GetEntry("Exit");
                  popup->HideEntry(exit->GetEntryId());
               }
            }
            ShowMenu(menu);
            return;
         }
      }
   }
   if (fActMenuBar != fMenuBar)
      ShowMenu(fMenuBar);
}

//______________________________________________________________________________
void TGNewBrowser::RecursiveReparent(TGPopupMenu *popup)
{
   TGMenuEntry *entry = 0;
   TIter next(popup->GetListOfEntries());
   while ((entry = (TGMenuEntry *)next())) {
      if (entry->GetPopup()) {
         RecursiveReparent(entry->GetPopup());
      }
   }
   popup->ReparentWindow(gClient->GetDefaultRoot());
}

//______________________________________________________________________________
/*
void TGNewBrowser::ReallyDelete()
{
   // Really delete the browser and the this GUI.

   gInterpreter->DeleteGlobal(fBrowser);
   delete fBrowser;    // will in turn delete this object
}
*/

//______________________________________________________________________________
TGFileBrowser* TGNewBrowser::MakeFileBrowser()
{
   // Create a file-browser. Caller should provide
   // Start/StopEmbedding() calls and populate the new browser.

   TBrowserImp    imp;
   TBrowser      *tb = new TBrowser("Pipi", "Strel", &imp);
   TGFileBrowser *fb = new TGFileBrowser(tb, "File browser", 200, 500);
   tb->SetBrowserImp((TBrowserImp *)fb);
   fb->SetBrowser(tb);
   fb->SetNewBrowser(this);
   return fb;
}

//______________________________________________________________________________
void TGNewBrowser::InitPlugins()
{
   // Initialize default plugins. Could be also of the form:
   // StartEmbedding(0);
   // TPluginHandler *ph;
   // ph = gROOT->GetPluginManager()->FindHandler("TGClassBrowser");
   // if (ph && ph->LoadPlugin() != -1) {
   //    ph->ExecPlugin(3, gClient->GetRoot(), 200, 500);
   // }
   // StopEmbedding();

   // --- Left vertical area

   // File browser plugin...
   StartEmbedding(0);
   //gROOT->ProcessLine(Form("new TGFileBrowser((const TGWindow *)0x%lx, 200, 500)", 
   //                   gClient->GetRoot()));
   {
      TGFileBrowser *fb = MakeFileBrowser();
      fb->BrowseObj(gROOT);
      fb->AddFSDirectory("/");
      fb->GotoDir(gSystem->WorkingDirectory());
      fb->Show();
   }
   StopEmbedding();
   SetTabTitle("Files", 0);

   // Class browser plugin
   StartEmbedding(0);
   gROOT->ProcessLine(Form("new TGClassBrowser((const TGWindow *)0x%lx, 200, 500)", 
                      gClient->GetRoot()));
   StopEmbedding();
   SetTabTitle("Classes", 0);

   // --- Right main area

   // Editor plugin...
   StartEmbedding(1);
   gROOT->ProcessLine(Form("new TGTextEditor((const char *)0, (const TGWindow *)0x%lx)", 
                      gClient->GetRoot()));
   StopEmbedding();
   SetTabTitle("Editor 1", 1);

   // HTML plugin...
   gSystem->Load("libGuiHtml");
   if (gSystem->Load("libRHtml") >= 0) {
      StartEmbedding(1);
      gROOT->ProcessLine(Form("new TGHtmlBrowser(\"http://root.cern.ch/root/html/ClassIndex.html\", \
                              (const TGWindow *)0x%lx)", gClient->GetRoot()));
      StopEmbedding();
      SetTabTitle("HTML", 1);
   }

   // Canvas plugin...
   StartEmbedding(1);
   gROOT->ProcessLine("new TCanvas(\"BrowserCanvas1\", \"Browser Canvas 1\")");
   StopEmbedding();
   SetTabTitle("Canvas 1", 1);

#if 0
   // GLViewer plugin...
   StartEmbedding(1);
   gROOT->ProcessLine(Form("new TGLSAViewer((TGFrame *)0x%lx, 0);", gClient->GetRoot()));
   StopEmbedding();
   SetTabTitle("GL Viewer", 1);

   // PROOF plugin...
   StartEmbedding(1);
   gROOT->ProcessLine("new TSessionViewer();");
   StopEmbedding();
   SetTabTitle("PROOF", 1);
#endif

   // --- Right bottom area

   // Command plugin...
   StartEmbedding(2);
   gROOT->ProcessLine(Form("new TGCommandPlugin((const TGWindow *)0x%lx, 700, 300)", 
                      gClient->GetRoot()));
   StopEmbedding();
   SetTabTitle("Command", 2);

   // --- Select first tab everywhere
   SetTab(0, 0);
   SetTab(1, 0);
   SetTab(2, 0);
}

/*
//______________________________________________________________________________
TBrowserImp *TGNewBrowser::NewBrowser(TBrowser *b, const char *title, UInt_t width, UInt_t height)
{
   TGNewBrowser *browser = new TGNewBrowser(b, title, width, height);
   return (TBrowserImp *)browser;
}

//______________________________________________________________________________
TBrowserImp *TGNewBrowser::NewBrowser(TBrowser *b, const char *title, Int_t x, Int_t y, UInt_t width, UInt_t height)
{
   TGNewBrowser *browser = new TGNewBrowser(b, title, x, y, width, height);
   return (TBrowserImp *)browser;
}
*/
