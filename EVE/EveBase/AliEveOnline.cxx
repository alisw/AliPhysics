//
//  AliEveOnline.cxx
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#include "AliEveOnline.h"
#include "AliEveGeomGentle.h"
#include "AliEveEventManager.h"
#include "AliEveEventManagerEditor.h"
#include "AliEveMultiView.h"
#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"

#include <TROOT.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TEveMacro.h>
#include <TGTab.h>
#include <TGFileBrowser.h>
#include <TInterpreter.h>
#include <TPRegexp.h>
#include <TFolder.h>
#include <TSystemDirectory.h>

#include <iostream>

using namespace std;

AliEveOnline::AliEveOnline()
{
    cout<<"Creating AliEveOnline"<<endl;
    
    // set OCDB path:
    AliEveEventManager::SetCdbUri("local:///local/cdb");         // current OCDB snapshot
    
    cout<<"Creating multiview...";
    AliEveMultiView *multiView = new AliEveMultiView(kTRUE);
    cout<<"created"<<endl;
    
    Info("alieve_init", "Adding standard macros.");
    TString  hack = gSystem->pwd(); // Problem with TGFileBrowser cding
    InitImportMacros();
    gSystem->cd(hack);
    cout<<"Standard macros added"<<endl;
    
    new AliEveEventManager("online", -1);
    gEve->AddEvent(AliEveEventManager::GetMaster());
    cout<<"Event manager created"<<endl;
    
    TEveUtil::AssertMacro("VizDB_scan.C");
    gSystem->ProcessEvents();
    cout<<"VizDB_scan loaded"<<endl;
    TEveBrowser *browser = gEve->GetBrowser();
    browser->ShowCloseTab(kFALSE);
    cout<<"browser created"<<endl;
    
    AliEveGeomGentle *geomGentle = new AliEveGeomGentle();
    
    multiView->InitGeomGentle(geomGentle->GetGeomGentle(),
                              geomGentle->GetGeomGentleRphi(),
                              geomGentle->GetGeomGentleRhoz(),
                              geomGentle->GetGeomGentleRhoz());
    
    multiView->InitGeomGentleTrd(geomGentle->GetGeomGentleTRD());
    multiView->InitGeomGentleMuon(geomGentle->GetGeomGentleMUON(), kFALSE, kFALSE, kTRUE);
    
    printf("============ Setting macro executor ============\n");
    
    AliEveMacroExecutor *exec = AliEveEventManager::GetMaster()->GetExecutor();
    printf("exec created\n");
    
    // default appearance:
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category",  "esd_tracks.C", "esd_tracks_by_category",  "", kTRUE));
    
    // preset for cosmics:
    //exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category",  "esd_tracks.C", "esd_tracks_by_category",  "kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kFALSE", kTRUE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD AD"   , "ad_esd.C", "ad_esd", "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD EMCAL", "emcal_esdcells.C", "emcal_esdcells", "", kTRUE));

    cout<<"macros added to exec"<<endl;
    
    //============================================================================
    // Final GUI setup
    //============================================================================
    
    browser->GetTabRight()->SetTab(1);
    browser->StartEmbedding(TRootBrowser::kBottom);
    new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
    browser->StopEmbedding("EventCtrl");
    
    browser->MoveResize(0, 0, gClient->GetDisplayWidth(),gClient->GetDisplayHeight() - 32);
    
    gEve->FullRedraw3D(kTRUE);
    gSystem->ProcessEvents();
    
    // move and rotate sub-views
    TGLViewer *glv1 = multiView->Get3DView()->GetGLViewer();
    TGLViewer *glv2 = multiView->GetRPhiView()->GetGLViewer();
    TGLViewer *glv3 = multiView->GetRhoZView()->GetGLViewer();
    TGLViewer *glv4 = multiView->GetMuonView()->GetGLViewer();
    
    glv1->CurrentCamera().RotateRad(-0.4, 0.6);
    glv2->CurrentCamera().Dolly(90, kFALSE, kFALSE);
    glv3->CurrentCamera().Dolly(2300, kFALSE, kFALSE);
    glv4->CurrentCamera().Dolly(1, kFALSE, kFALSE);
    
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();
    gEve->Redraw3D(kTRUE);
    
    // set autoload by default
    AliEveEventManager::GetMaster()->SetAutoLoad(true);
    AliEveEventManager::GetMaster()->SetSaveViews(true);

}

AliEveOnline::~AliEveOnline()
{
    
}

void AliEveOnline::InitImportMacros()
{
    // Put macros in the list of browsables, add a macro browser to
    // top-level GUI.
    
    TString macdir("$(ALICE_ROOT)/EVE/alice-macros");
    
    if (gSystem->Getenv("ALICE_ROOT") != 0)
    {
        gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
        gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
        gSystem->ExpandPathName(macdir);
    }
    
    
    TFolder* f = gEve->GetMacroFolder();
    void* dirhandle = gSystem->OpenDirectory(macdir.Data());
    if (dirhandle != 0)
    {
        char* filename;
        TPMERegexp re("\\.C$");
        TObjArray names;
        while ((filename = (char*)(gSystem->GetDirEntry(dirhandle))) != 0)
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