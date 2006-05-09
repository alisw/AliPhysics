// $Header$

void alieve_init(const Text_t* path=".", Int_t event=0,
		  Bool_t use_runloader=true, Bool_t use_esd=true)
{

  // Set-up environment, load libraries.

  Reve::SetupEnvironment();

  gROOT->SetMacroPath(Form("%s:%s/alice-macros:%s/macros",
			   gROOT->GetMacroPath(),
			   gSystem->Getenv("REVESYS"),
			   gSystem->Getenv("ALICE_ROOT")));
  gInterpreter->AddIncludePath(Form("%s/macros", gSystem->Getenv("ALICE_ROOT")));

  Reve::AssertMacro("alieve_loadlibs.C");
  gSystem->Load("libAlieve.so");


  // Put macros in the list of browsables, spawn a browser.

  TFolder* f = new TFolder("ALICE EVE", "Visualization macros");
  TString macdir("$(REVESYS)/alice-macros");
  gSystem->ExpandPathName(macdir);

  void* dirhandle = gSystem->OpenDirectory(macdir.Data());
  if(dirhandle != 0) {
    char* filename;
    TPRegexp re("\.C$");
    while((filename = gSystem->GetDirEntry(dirhandle)) != 0) {
      if(re.Match(filename)) {
	printf("Adding macro '%s'\n", filename);
	f->Add(new TMacro(Form("%s/%s", macdir.Data(), filename)));
      }
    }
  }
  gSystem->FreeDirectory(dirhandle);

  gROOT->GetListOfBrowsables()->Add(f);
  gROOT->GetListOfBrowsables()->Add
    (new TSystemDirectory(macdir.Data(), macdir.Data()));

  new TBrowser;

  // Open event

  Alieve::Event::Initialize(use_runloader, use_esd);

  printf("Opening event %d from '%s' ...", event, path); fflush(stdout);
  Alieve::gEvent = new Alieve::Event(path, event);
  printf(" done.\n");
  gReve->AddEvent(Alieve::gEvent);
}
