// $Header$

void alieve_init(const Text_t* path=".", Int_t event=0,
		 Bool_t use_runloader=kTRUE, Bool_t use_esd=kTRUE,
		 Bool_t avoid_exceptions_on_open=kTRUE)
{

  // Set-up environment, load libraries.

  Reve::SetupEnvironment();
  // alieve executable linked against ALICE nad EVE shared libraries.
  // Reve::AssertMacro("alieve_loadlibs.C");


  // Put macros in the list of browsables, spawn a browser.

  TString macdir("$(REVESYS)/alice-macros");
  gSystem->ExpandPathName(macdir);

  TFolder* f = gReve->GetMacroFolder();
  void* dirhandle = gSystem->OpenDirectory(macdir.Data());
  if(dirhandle != 0) {
    char* filename;
    TPRegexp re("\.C$");
    while((filename = gSystem->GetDirEntry(dirhandle)) != 0) {
      if(re.Match(filename)) {
	printf("Adding macro '%s'\n", filename);
	f->Add(new Reve::RMacro(Form("%s/%s", macdir.Data(), filename)));
      }
    }
  }
  gSystem->FreeDirectory(dirhandle);

  gROOT->GetListOfBrowsables()->Add
    // (new TSystemDirectory("alice-macros", macdir.Data())); // !!!! this spits blood, but then works
    (new TSystemDirectory(macdir.Data(), macdir.Data()));

  new TBrowser;

  Reve::AssertMacro("region_marker.C");

  // Open event
  if(path != 0) {
    Alieve::Event::Initialize(use_runloader, use_esd, avoid_exceptions_on_open);

    printf("Opening event %d from '%s' ...", event, path); fflush(stdout);
    Alieve::gEvent = new Alieve::Event(path, event);
    printf(" done.\n");
    gReve->AddEvent(Alieve::gEvent);
  }
}
