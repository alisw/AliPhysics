// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file alieve_init.C

#if !defined(__CINT__) || defined(__MAKECINT__)
//#include <AliQAHistViewer.h>

#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TInterpreter.h>
#include <TMath.h>
#include <TGListTree.h>
#include <TEveVSDStructs.h>
#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEnv.h>
#include <TEveWindowManager.h>
#include <TGTab.h>
#include <TTimeStamp.h>
#include <TGFileBrowser.h>
#include <TPRegexp.h>
#include <TGrid.h>
#include <TFolder.h>
#include <TEveMacro.h>
#include <TSystemDirectory.h>

#include <AliESDtrackCuts.h>
//#include <AliPWG0Helper.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliESDtrack.h>
#include <AliESDfriendTrack.h>
#include <AliExternalTrackParam.h>
#include <AliEveTrack.h>
#include <AliEveTrackCounter.h>
#include <AliEveMagField.h>
#include <AliEveEventManager.h>
#include <AliEveEventManagerEditor.h>
#include <AliEveMultiView.h>
#include <AliEveMacroExecutor.h>
#include <AliEveMacro.h>
#include <AliEveEventSelectorWindow.h>
#include <AliEveDataSourceOffline.h>
#include <AliEveTrackFitter.h>
#include <AliCDBManager.h>

#endif

void alieve_init_import_macros();

void alieve_init(const TString& cdburi = "",
                 const TString& path   = ".", Int_t event=0,Bool_t,
                 const Text_t* esdfile = 0,
                 const Text_t* aodfile = 0,
                 const Text_t* rawfile = 0,
                 Bool_t assert_runloader = kFALSE,
                 Bool_t assert_esd       = kFALSE,
                 Bool_t assert_aod       = kFALSE,
                 Bool_t assert_raw       = kFALSE)
{
    if (cdburi.IsNull() && ! AliCDBManager::Instance()->IsDefaultStorageSet())
    {
        gEnv->SetValue("Root.Stacktrace", "no");
        Fatal("alieve_init.C", "OCDB path MUST be specified as the first argument.");
    }
    
    AliEveEventManager::SetCdbUri(cdburi);
    
    Info("alieve_init", "Adding standard macros.");
    TString  hack = gSystem->pwd(); // Problem with TGFileBrowser cding
    alieve_init_import_macros();
    gSystem->cd(hack);
    
    TEveUtil::AssertMacro("VizDB_scan.C");
    
    gSystem->ProcessEvents();
    
    // Open event
    if (path.BeginsWith("alien:") || ! cdburi.BeginsWith("local:"))
    {
        if (gGrid != 0)
        {
            Info("alieve_init", "TGrid already initializied. Skiping checks and initialization.");
        }
        else
        {
            Info("alieve_init", "AliEn requested - connecting.");
            if (gSystem->Getenv("GSHELL_ROOT") == 0)
            {
                Error("alieve_init", "AliEn environment not initialized. Aborting.");
                gSystem->Exit(1);
            }
            if (TGrid::Connect("alien") == 0)
            {
                Error("alieve_init", "TGrid::Connect() failed. Aborting.");
                gSystem->Exit(1);
            }
        }
    }
    
    AliEveDataSourceOffline *dataSource = (AliEveDataSourceOffline*)AliEveEventManager::Instance()->GetDataSourceOffline();
    
    dataSource->SetFilesPath(path);
    
    Info("alieve_init", "Opening event %d from '%s' ...", event, path.Data());
    TString name("Event"); // CINT has trouble with direct "Event".
    //  new AliEveEventManager(name, event);
    gEve->AddEvent(AliEveEventManager::Instance());
}

void alieve_init_import_macros()
{
    // Put macros in the list of browsables, add a macro browser to
    // top-level GUI.
    
    TString macdir("$(ALICE_ROOT)/EVE/macros");
    gSystem->ExpandPathName(macdir);
    
    TFolder* f = gEve->GetMacroFolder();
    void* dirhandle = gSystem->OpenDirectory(macdir.Data());
    if (dirhandle != 0)
    {
        const char* filename;
        TPMERegexp re("\\.C$");
        TObjArray names;
        while ((filename = gSystem->GetDirEntry(dirhandle)) != 0)
        {
            if (re.Match(filename))
                names.AddLast(new TObjString(filename));
        }
        names.Sort();
        
        for (Int_t ii=0; ii<names.GetEntries(); ++ii)
        {
            TObjString * si = (TObjString*) names.At(ii);
            f->Add(new TEveMacro(Form("%s/%s", macdir.Data(), (si->GetString()).Data())));
        }
    }
    gSystem->FreeDirectory(dirhandle);
    
    gROOT->GetListOfBrowsables()->Add(new TSystemDirectory(macdir.Data(), macdir.Data()));
    
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
