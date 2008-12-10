// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void alieve_init(const Text_t* path   = ".", Int_t event=0,
                 const Text_t* esdfile = 0,
                 const Text_t* aodfile = 0,
                 const Text_t* rawfile = 0,
		 const Text_t* cdburi  = 0,
		 Bool_t assert_runloader = kFALSE,
                 Bool_t assert_esd       = kFALSE,
                 Bool_t assert_aod       = kFALSE,
                 Bool_t assert_raw       = kFALSE)
{
  Info("alieve_init", "Adding standard macros.");
  TString  hack = gSystem->pwd(); // Problem with TGFileBrowser cding
  alieve_init_import_macros();
  gSystem->cd(hack);

  alieve_init_basic_vizdb();
  // Temporarily assert also default vizdb.
  TEveUtil::AssertMacro("VizDB_scan.C");

  gSystem->ProcessEvents();

  AliEveEventManager::SetESDFileName(esdfile);
  AliEveEventManager::SetRawFileName(rawfile);
  AliEveEventManager::SetCdbUri(cdburi);
  AliEveEventManager::SetAssertElements(assert_runloader, assert_esd,
					assert_aod, assert_raw);

  // Open event
  if (path != 0)
  {
    Info("alieve_init", "Opening event %d from '%s' ...", event, path);
    TString name("Event"); // CINT has trouble with direct "Event".
    new AliEveEventManager(name, path, event);
    gEve->AddEvent(AliEveEventManager::GetMaster());
  }
}

void alieve_init_import_macros()
{
  // Put macros in the list of browsables, add a macro browser to
  // top-level GUI.

  TString macdir("$(ALICE_ROOT)/EVE/alice-macros");
  gSystem->ExpandPathName(macdir);

  TFolder* f = gEve->GetMacroFolder();
  void* dirhandle = gSystem->OpenDirectory(macdir.Data());
  if (dirhandle != 0)
  {
    char* filename;
    TPMERegexp re("\\.C$");
    std::list<string> names; // This form understood by cint (fails with std::string).
    while ((filename = gSystem->GetDirEntry(dirhandle)) != 0)
    {
      if (re.Match(filename))
	names.push_back(filename);
    }
    names.sort();
    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    // f->Add(new TEveMacro(Form("%s/%s", macdir.Data(), filename)));
    char fullName[1000];
    for (std::list<string>::iterator si=names.begin(); si!=names.end(); ++si)
    {
      sprintf(fullName,"%s/%s", macdir.Data(), si->c_str());
      f->Add(new TEveMacro(fullName));
    }
  }
  gSystem->FreeDirectory(dirhandle);

  gROOT->GetListOfBrowsables()->Add
    // (new TSystemDirectory("alice-macros", macdir.Data())); // !!!! this spits blood, but then works
    (new TSystemDirectory(macdir.Data(), macdir.Data()));

  {
    TEveBrowser   *br = gEve->GetBrowser();
    TGFileBrowser *fb = 0;
    fb = br->GetFileBrowser();
    fb->GotoDir(macdir);
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
}

void alieve_init_basic_vizdb()
{
  TEvePointSet* ps;

  ps = new TEvePointSet();
  ps->SetMarkerColor(4);
  ps->SetMarkerSize(0.2);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("Clusters", ps);
}
