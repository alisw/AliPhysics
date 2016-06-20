// $Id: alieve_init.C 30728 2009-01-22 18:14:34Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TFolder.h"
#include "TPRegexp.h"
#include "TEveUtil.h"

#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveBrowser.h"

#include "TSystemDirectory.h"
#include "TGFileBrowser.h"
#include "TEveMacro.h"
#endif

using namespace std;

typedef list<string> StringList;


void alieve_init_import_macros()
{
  // Put macros in the list of browsables, add a macro browser to
  // top-level GUI.
  typedef list<string> StringList;

  TString macdir("$(ALICE_ROOT)/EVE/macros");
  gSystem->ExpandPathName(macdir);

  TFolder* f = gEve->GetMacroFolder();
  void* dirhandle = gSystem->OpenDirectory(macdir.Data());
  if (dirhandle != 0)
  {
    char* filename;
    TPMERegexp re("\\.C$");
    StringList names; // This form understood by cint (fails with std::string).
    while ((filename = (char*)gSystem->GetDirEntry(dirhandle)) != 0)
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
    for (StringList::iterator si=names.begin(); si!=names.end(); ++si)
    {
      sprintf(fullName,"%s/%s", macdir.Data(), si->c_str());
      f->Add(new TEveMacro(fullName));
    }
  }
  gSystem->FreeDirectory(dirhandle);

  gROOT->GetListOfBrowsables()->Add
    // (new TSystemDirectory("macros", macdir.Data())); // !!!! this spits blood, but then works
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

void hlt_alieve_init()
{
  typedef list<string> StringList;

  Info("alieve_init", "Adding standard macros.");
  TString  hack = gSystem->pwd(); // Problem with TGFileBrowser cding
  alieve_init_import_macros();
  gSystem->cd(hack);

  alieve_init_basic_vizdb();
  // Temporarily assert also default vizdb.
  TEveUtil::AssertMacro("VizDB_scan.C");

  gSystem->ProcessEvents();
}

