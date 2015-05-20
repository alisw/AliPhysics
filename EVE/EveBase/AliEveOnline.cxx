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

AliEveOnline::AliEveOnline(bool storageManager)
{
    cout<<"Creating AliEveOnline"<<endl;
 
    //-----------------------------------------------------------------------------------------
    //  Set all preferences here
    //
    Color_t colorTRD = kGray;      // color of TRD modules
    Color_t colorMUON = kGray;     // color of MUON modules
    
    bool customPreset = true;            // should one of the following custom presets be used
    Color_t colors[9] = {kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan};
    Width_t widths[9] = {3,3,3,3,3,3,3,3,3};
    bool dashBad = true;
    
//    Color_t colors[9] = {kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen}; // preset for cosmics

    bool saveViews = true;          // should screenshot be saved and sent to ALICE LIVE
    
    //
    //-----------------------------------------------------------------------------------------
    
    // set OCDB path:
    AliEveEventManager::SetCdbUri("local:///local/cdb");         // current OCDB snapshot
    
    cout<<"Creating multiview...";
    AliEveMultiView *multiView = new AliEveMultiView(kTRUE);
    cout<<"created"<<endl;
    
    cout<<"Adding standard macros...";
    InitImportMacros();
    cout<<"added"<<endl;
    
    cout<<"Creating event manager...";
    new AliEveEventManager("online", -1);
    AliEveEventManager *man = AliEveEventManager::GetMaster();
    gEve->AddEvent(man);
    cout<<"created"<<endl;
    
    TEveUtil::AssertMacro("VizDB_scan.C");
    gSystem->ProcessEvents();
    cout<<"VizDB_scan loaded"<<endl;
    TEveBrowser *browser = gEve->GetBrowser();
    browser->ShowCloseTab(kFALSE);
    cout<<"Browser created"<<endl;
    
    AliEveGeomGentle *geomGentle = new AliEveGeomGentle();
    
    multiView->InitGeomGentle(geomGentle->GetGeomGentle(),
                              geomGentle->GetGeomGentleRphi(),
                              geomGentle->GetGeomGentleRhoz(),
                              geomGentle->GetGeomGentleRhoz());
    
    multiView->InitGeomGentleTrd(geomGentle->GetGeomGentleTRD(colorTRD));
    multiView->InitGeomGentleMuon(geomGentle->GetGeomGentleMUON(true,colorMUON), kFALSE, kFALSE, kTRUE);
    
    cout<<"============ Setting macro executor ============\n"<<endl;;
    AliEveMacroExecutor *exec = AliEveEventManager::GetMaster()->GetExecutor();
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD AD"   , "ad_esd.C", "ad_esd", "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD EMCAL", "emcal_esdcells.C", "emcal_esdcells", "", kTRUE));
    cout<<"macros added to exec"<<endl;
    
    //============================================================================
    // Final GUI setup
    //============================================================================
    
    browser->GetTabRight()->SetTab(1);
    browser->StartEmbedding(TRootBrowser::kBottom);
    new AliEveEventManagerWindow(man,storageManager);
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

    if(customPreset)
    {
        man->SetESDcolors(colors);
        man->SetESDwidths(widths);
        man->SetESDdashBad(dashBad);
    }

    man->SetSaveViews(saveViews);
    man->SetAutoLoad(true);     // set autoload by default
}

AliEveOnline::~AliEveOnline()
{
}

void AliEveOnline::InitImportMacros()
{
    // Put macros in the list of browsables, add a macro browser to
    // top-level GUI.
    
    TString  hack = gSystem->pwd(); // Problem with TGFileBrowser cding
    
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
    gSystem->cd(hack);
}