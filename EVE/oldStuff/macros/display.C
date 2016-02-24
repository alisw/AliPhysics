#ifndef __CINT__
# include <TGFileDialog.h>
# include <TGFileBrowser.h>
# include <TControlBar.h>
# include <TString.h>
# include <TSystem.h>
# include <TROOT.h>
# include <TEveManager.h>
# include <TSystemDirectory.h>
# include <TEveMacro.h>
# include <TEveBrowser.h>
# include <TFolder.h>
# include <TList.h>
# include <TApplication.h>
# include "aliroot/EVE/EveBase/AliEveEventManager.h"
# include <iostream>
#endif

class Display
{
public:
  //____________________________________________________________________
  /** Constructor 
      @param file File or directory to read data from */
  Display(const char* file="./")
    : fCb(0), 
      fEnableHits(kFALSE),
      fEnableDigits(kFALSE),
      fEnableRaw(kFALSE),
      fEnableESD(kFALSE),
      fHasHits(kFALSE),
      fHasDigits(kFALSE),
      fHasRaw(kFALSE),
      fHasESD(kFALSE)
  {
    AliEveEventManager::SetCdbUri(TString("local://$ALICE_ROOT/OCDB"));
    LoadMacros();
    Setup(file, false);
    SetupSelect();
    SetupControl();
  }
  //____________________________________________________________________
  /** Setup the control bar */
  void SetupControl()
  {
    fCb = new TControlBar("vertical", "EVE Controls");
    const char* buts[] = { "Load ...", 
			   "Select ...",
			   "Next event", 
			   "Previous event", 
			   "Hits", 
			   "Digits", 
			   "Raw", 
			   "ESD", 
			   "Tracks", 
			   "Quit", 
			   0 }; 
    const char* meths[] = { "LoadFile", 
			    "Select",
			    "NextEvent", 
			    "PrevEvent", 
			    "ToggleHits", 
			    "ToggleDigits", 
			    "ToggleRaw", 
			    "ToggleESD", 
			    "ToggleTracks", 
			    "Quit", 
			    0 };
    const char* hints[] = { "Specify input" , 
			    "Select which detectors to show",
			    "Get next event", 
			    "Get previous event", 
			    "Toggle display of hits", 
			    "Toggle display of simulated digits", 
			    "Toggle display of raw digits", 
			    "Toggle display of event summary data", 
			    "Toggle display of tracks", 
			    "Quit", 
			    0 };
			    
    char** but  = const_cast<char**>(buts);
    char** meth = const_cast<char**>(meths);
    char** hint = const_cast<char**>(hints);
    while (*but) { 
      fCb->AddButton(*but, Form("Display::Instance()->%s()", *meth), *hint);
      but++;
      meth++;
      hint++;
    }
    fCb->Show();
  }
  //____________________________________________________________________
  /** Setup the selection bar */
  void SetupSelect()
  {
    fSelect = new TControlBar("vertical", "Detector selection");
    char* dets[] = { "emcal", 
		     "fmd", 
		     "its", 
		     "pmd", 
		     "t0", 
		     "tof", 
		     "tpc", 
		     "trd", 
		     "vzero",
		     "kine", 
		     "esd",
		     0 };
    char** det = dets;
    while (*det) { 
      // TString cmd(Form("Display::Instance()->Enable(\"%s\")", *det));
      // TString help(Form("Enable display of %s data", name.Data()));
      fSelect->AddButton(*det, 
			 Form("Display::Instance()->Enable(\"%s\")", *det),
			 Form("Enable display of %s data", *det));
      det++;
    }
    std::cout << "Adding \"Hide\" to select" << std::endl;
    fSelect->AddButton("Hide", "Display::Instance()->HideSelect()", 
		       "Hide the detector selection menu");
  }
  //____________________________________________________________________
  /** Load macros used by EVE */
  void LoadMacros()
  {
    std::cout << "Loading macros ... " << std::flush;
    TString savdir(gSystem->pwd());
    TString macdir("$(ALICE_ROOT)/EVE/macros");
    gSystem->ExpandPathName(macdir);
  
    TFolder*          f     = gEve->GetMacroFolder();
    TSystemDirectory* dir   = new TSystemDirectory(macdir.Data(), 
						   macdir.Data());
    TList*            files = dir->GetListOfFiles();
    files->Sort();
    TIter             next(files);
    TSystemFile*      file  = 0;
    while ((file = static_cast<TSystemFile*>(next()))) { 
      if (file->IsDirectory()) continue;
      TString name(gSystem->ConcatFileName(macdir.Data(),file->GetName()));
      if (!name.EndsWith(".C")) continue;
      f->Add(new TEveMacro(name.Data()));
    }
    
    gROOT->GetListOfBrowsables()->Add(dir);
  
    {
      TEveBrowser*   b  = gEve->GetBrowser();
      TGFileBrowser* fb = b->GetFileBrowser();
      fb->GotoDir(macdir);
      { 
	b->StartEmbedding(0);
	fb = b->MakeFileBrowser();
	fb->BrowseObj(f);
	fb->Show();
	b->StopEmbedding(0);
	b->SetTabTitle("Macros", 0);
	b->SetTab(0,0);
      }
    }
    gSystem->cd(savdir);
    std::cout << "done" << std::endl;
  }
  //____________________________________________________________________
  /** Set-up load location, etc. 
      @param file File or directory to read data from 
      @param refresh Force a refresh */
  void Setup(const char* file, bool refresh)
  {
    TString fileName(gSystem->BaseName(file));
    TString dirName(gSystem->DirName(file));
    
    if (fileName.Contains("ESD")) { 
      std::cout << "Adding ESD file " << fileName << std::endl;
      AliEveEventManager::Instance()->GetDataSourceOffline()->SetESDFileName(fileName);
    }
    else if (fileName.EndsWith(".root") || fileName.EndsWith(".raw")) { 
      std::cout << "Adding raw file " << fileName << std::endl;
      AliEveEventManager::SetRawFileName(fileName);
    }
    else if (fileName.IsNull()) { 
      std::cout << "No file given!" << std::endl;
      // std::cout << "Adding raw directory " << dirName.Data() << std::endl;
      // AliEveEventManager::SetRawFileName(dirName.Data());
    }
    else {
      std::cerr << "Don't know how to deal with '" << fileName << "'" 
		<< std::endl;
      return;
    }
    std::cout << "Opening " << fileName << " (" << dirName << ")" << std::endl;
    
    if (AliEveEventManager::Instance()) delete AliEveEventManager::Instance();
    TString eventName("Event"); // CINT has trouble with direct "Event".
    /* AliEveEventManager::Instance() =*/
    new AliEveEventManager(eventName, dirName, 0);
    gEve->AddEvent(AliEveEventManager::Instance());
    
    if (refresh) Refresh();
  }
  //____________________________________________________________________
  /** Executed in response to hitting the "Load ..." button. */
  void 
  LoadFile()
  {
    fCb->SetButtonState("Load ...", 1);
    TGFileInfo fi;
    fi.fIniDir = StrDup(gSystem->pwd());
    new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);
    Setup(fi.fFilename, true);
    fCb->SetButtonState("Load ...", 0);
  }
  //____________________________________________________________________
  /** Pop-up or down the selection menu */
  void 
  Select()
  {
    fSelect->Show();
    fCb->SetButtonState("Select ...", 1);
  }
  //____________________________________________________________________
  /** Enable a detector or type of display */
  void 
  Enable(const char* what)
  {
    TObject* o = fEnabled.FindObject(what);
    //TString n(what);
    // n.ToUpper();
    std::cout << "Enable " << what << "?" << std::endl;
    if (o) { 
      fEnabled.Remove(o);
      fSelect->SetButtonState(what, 0);
      delete o;
    }
    else { 
      fEnabled.Add(new TObjString(what));
      fSelect->SetButtonState(what, 1);
    }
    Refresh();
  }
  //____________________________________________________________________
  /** Pop-down the selection menu */
  void 
  HideSelect()
  {
    fCb->SetButtonState("Select ...", 0);
    fSelect->Hide();
  }
  //____________________________________________________________________
  /** Prepare for a new event */
  void 
  Prepare()
  {
    fHasHits   = kFALSE;
    fHasDigits = kFALSE;
    fHasRaw    = kFALSE;
    fHasESD    = kFALSE;
    fHasTracks = kFALSE;
  }
  //____________________________________________________________________
  /** Get the next event */
  void 
  NextEvent()
  {
    Prepare();
    std::cout << "Getting next event, please wait ... " << std::flush;
    gROOT->Macro("event_next.C");
    std::cout << "done" << std::endl;
    Refresh();
  }
  //____________________________________________________________________
  /** Get the previous event */
  void 
  PrevEvent()
  {
    Prepare();
    std::cout << "Getting previous event, please wait ... " << std::flush;
    gROOT->Macro("event_prev.C"); 
    std::cout << "done" << std::endl;
    Refresh();
  }
  //____________________________________________________________________
  /** Reload the current event */
  void
  Reload()
  {
    Prepare();
    Int_t event = AliEveEventManager::Instance()->GetEventId();
    std::cout << "Getting event " << event 
	      << ", please wait ... " << std::flush;
    gROOT->Macro(Form("event_goto.C(%d)", event));
    std::cout << "done" << std::endl;
    Refresh();
  }
  //__________________________________________________________________
  /** Refresh the display of something
      @param enabled Whether this something is enabled or not 
      @param has     Whether we have already shown this data 
      @param what    What data to show */
  void 
  RefreshOne(Bool_t enabled, Bool_t& has, const char* what)
  {
    if (!enabled || has) return;
    std::cout << what << "..." << std::flush;
    TIter next(&fEnabled);
    TObject* o = 0;
    while ((o = next())) { 
      TMacro* macro = gEve->GetMacro(Form("%s_%s", o->GetName(), what));
      if (macro) macro->Exec();
    }
    has = kTRUE;
  }
  
  //____________________________________________________________________
  /** Refresh the display of the data */
  void 
  Refresh()
  {
    std::cout << "Drawing ... " << std::flush;
    RefreshOne(fEnableTracks, fHasTracks, "tracks");
    RefreshOne(fEnableHits,   fHasHits,   "hits");
    RefreshOne(fEnableDigits, fHasDigits, "digits");
    RefreshOne(fEnableRaw,    fHasRaw,    "raw");
    RefreshOne(fEnableESD,    fHasESD,    "esd");
    std::cout << "done" << std::endl;
  }
  //____________________________________________________________________
  /** Toggle the diplay of some data */
  void 
  Toggle(Bool_t& what, const char* name)
  {
    what = !what;
    fCb->SetButtonState(name, what ? 1 : 0);
    if (!what) Reload();
    Refresh();
  }
  //____________________________________________________________________
  /** Toggle the diplay the simulated hits */
  void ToggleHits() { Toggle(fEnableHits, "Hits"); } 
  //____________________________________________________________________
  /** Toggle the diplay the simulated digits */
  void ToggleDigits() { Toggle(fEnableDigits, "Digits"); }
  //____________________________________________________________________
  /** Toggle the diplay the raw digits */
  void ToggleRaw() { Toggle(fEnableRaw, "Raw"); }
  //____________________________________________________________________
  /** Toggle the diplay the Event Summary Data */
  void ToggleESD() { Toggle(fEnableESD, "ESD"); }
  //____________________________________________________________________
  /** Toggle the diplay of tracks */
  void ToggleTracks() { Toggle(fEnableTracks, "Tracks"); }
  //____________________________________________________________________
  /** Pop-down both menues */
  void Quit() 
  { 
    fCb->Hide();
    fSelect->Hide();
    fgInstance = 0;
    delete this;
  }
  //____________________________________________________________________
  /** Get the static instance */
  static Display* Instance(const char* file=0) { 
    if (!fgInstance) fgInstance = new Display(file);
    return fgInstance;
  }
protected:
  TControlBar* fCb;       // Control bar
  TControlBar* fSelect;   // Selection bar
  Bool_t fEnableHits;     // Whether simulated hits are enabled 
  Bool_t fEnableDigits;   // Whether simulated digits are enabled 
  Bool_t fEnableRaw;      // Whehter raw data digits are enabled
  Bool_t fEnableESD;      // Whehter Event Summary Data is enabled
  Bool_t fEnableTracks;   // Whether tracks are shown
  Bool_t fHasHits;        // Whether simulated hits are rendered
  Bool_t fHasDigits;      // Whether simulated digits are rendered
  Bool_t fHasRaw;         // Whether raw data digits are rendered
  Bool_t fHasESD;         // Whether Event Summary Data is rendered
  Bool_t fHasTracks;      // Whether tracks are rendered
  TList  fEnabled;
  static Display* fgInstance; // Static singleton instance 
};
Display* Display::fgInstance = 0;


//____________________________________________________________________
Display* display(const char* file="")
{
  // TGeoManager::Import(geom);
  // AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  // AliCDBManager::Instance()->SetRun(0);
  // AliGeomManager::LoadGeometry("geometry.root");

  Display* d = Display::Instance(file);
  return d;

}
//____________________________________________________________________
//
// EOF
//
