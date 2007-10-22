// $Header$

#ifndef __CINT_
#include <list>
#include <string>
#endif

void alieve_init(const Text_t* path   = ".", Int_t event=0,
		 const Text_t* cdburi = 0,
		 Bool_t assert_runloader=kFALSE, Bool_t assert_esd=kFALSE)
{
  using namespace std;

  // Set-up environment, load libraries.

  Reve::SetupEnvironment();

  // Put macros in the list of browsables, spawn a browser.

  Info("alieve_init", "Adding standard macros.");

  TString macdir("$(REVESYS)/alice-macros");
  gSystem->ExpandPathName(macdir);

  TFolder* f = gReve->GetMacroFolder();
  void* dirhandle = gSystem->OpenDirectory(macdir.Data());
  if(dirhandle != 0) {
    char* filename;
    TPRegexp re("\.C$");
    list<string> names;
    while((filename = gSystem->GetDirEntry(dirhandle)) != 0) {
      if(re.Match(filename)) {
	names.push_back(filename);
      }
    }
    names.sort();
    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    // f->Add(new Reve::RMacro(Form("%s/%s", macdir.Data(), filename)));
    char fullName[1000];
    for (list<string>::iterator si=names.begin(); si!=names.end(); ++si)
    {
      sprintf(fullName,"%s/%s", macdir.Data(), si->c_str());
      f->Add(new Reve::RMacro(fullName));
    }
  }
  gSystem->FreeDirectory(dirhandle);

  gROOT->GetListOfBrowsables()->Add
    // (new TSystemDirectory("alice-macros", macdir.Data())); // !!!! this spits blood, but then works
    (new TSystemDirectory(macdir.Data(), macdir.Data()));

  {
    Reve::RGBrowser *br = gReve->GetBrowser();
    TGFileBrowser   *fb = 0;
    fb = br->GetFileBrowser();
    fb->GotoDir("/alice-macros"); //macdir);
    {
      br->StartEmbedding(0);
      fb = br->MakeFileBrowser();
      fb->BrowseObj(f);
      fb->Show();
      br->StopEmbedding();
      br->SetTabTitle("Macros", 0);
      br->SetTab(0, 0);
    }
  }

  // Reve::AssertMacro("region_marker.C");
  
  gSystem->ProcessEvents();

  // Open event
  if(path != 0) {
    Alieve::Event::SetCdbUri(cdburi);
    Alieve::Event::SetAssertElements(assert_runloader, assert_esd);
    printf("Opening event %d from '%s' ...", event, path); fflush(stdout);
    Alieve::gEvent = new Alieve::Event(path, event);
    printf(" done.\n");
    gReve->AddEvent(Alieve::gEvent);
  }
}
