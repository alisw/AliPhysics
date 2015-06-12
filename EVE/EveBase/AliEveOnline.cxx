//
//  AliEveOnline.cxx
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#include "AliEveOnline.h"
#include "AliEveGeomGentle.h"
#include "AliEveEventManagerEditor.h"
#include "AliEveMultiView.h"
#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"

#include <TROOT.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TBrowser.h>
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
    
    Color_t colors[9] = {kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan};
    Width_t width = 3;
    bool dashNoRefit = true;
    bool drawNoRefit = true;
    
//    Color_t colors[9] = {kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen}; // preset for cosmics

    bool saveViews = true;          // should screenshot be saved and sent to ALICE LIVE
    
    //
    //-----------------------------------------------------------------------------------------
 
    
    cout<<"Creating event manager...";
    AliEveEventManager *man = new AliEveEventManager(AliEveEventManager::kSourceOnline);
    gEve->AddEvent(man);
    cout<<"created"<<endl;
    
    // set OCDB path:
    AliEveEventManager::SetCdbUri("local:///local/cdb");         // current OCDB snapshot
    
    cout<<"Creating multiview...";
    AliEveMultiView *multiView = new AliEveMultiView(false);
    cout<<"created"<<endl;
    
    cout<<"Adding standard macros...";
    InitImportMacros();
    cout<<"added"<<endl;
    
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
                              0/*geomGentle->GetGeomGentleRhoz()*/);
    
    multiView->InitGeomGentleTrd(geomGentle->GetGeomGentleTRD(colorTRD));
    multiView->InitGeomGentleMuon(geomGentle->GetGeomGentleMUON(true,colorMUON), kFALSE, kTRUE, kFALSE);
    
    cout<<"============ Setting macro executor ============\n"<<endl;;
    AliEveMacroExecutor *exec = AliEveEventManager::GetMaster()->GetExecutor();
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM TrackRef MUON", "muon_trackRefs.C", "muon_trackRefs", "kTRUE", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks MUON", "esd_muon_tracks.C", "esd_muon_tracks", "kTRUE,kFALSE", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD AD"   , "ad_esd.C", "ad_esd", "", kTRUE));
//    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "PT HISTO"   , "histo2d.C", "histo2d", "", kTRUE));
    
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD EMCAL", "emcal_esdclustercells.C", "emcal_esdclustercells", "", kTRUE));
    
    cout<<"macros added to exec"<<endl;
    
    //============================================================================
    // Final GUI setup
    //============================================================================
    
    browser->GetTabRight()->SetTab(1);
    browser->StartEmbedding(TRootBrowser::kBottom);
    new AliEveEventManagerWindow(man,storageManager);
    browser->StopEmbedding("EventCtrl");
    
//    browser->MoveResize(0, 0, gClient->GetDisplayWidth(),gClient->GetDisplayHeight() - 32);
    
    gEve->FullRedraw3D(kTRUE);
    gSystem->ProcessEvents();
    
    // move and rotate sub-views
    TGLViewer *glv1 = multiView->Get3DView()->GetGLViewer();
    TGLViewer *glv2 = multiView->GetRPhiView()->GetGLViewer();
    TGLViewer *glv3 = multiView->GetRhoZView()->GetGLViewer();
//    TGLViewer *glv4 = multiView->GetMuonView()->GetGLViewer();
    
    glv1->CurrentCamera().RotateRad(-0.4, 0.6);
    glv2->CurrentCamera().Dolly(1, kFALSE, kFALSE);
    glv3->CurrentCamera().Dolly(1, kFALSE, kFALSE);
//    glv4->CurrentCamera().Dolly(1, kFALSE, kFALSE);
    
    
//    TGFrame *f1 = multiView->Get3DView()->GetGUIFrame();
//    TGFrame *f2 = multiView->GetRPhiView()->GetGUIFrame();
//    TGFrame *f3 = multiView->GetRhoZView()->GetGUIFrame();
//    TGFrame *f4 = multiView->GetMuonView()->GetGUIFrame();
//    
//    int width = f1->GetWidth()+f2->GetWidth();
//    int heigth = f1->GetHeight();
//    
////    TGWindow *superParent = (TGWindow*)f1->GetParent()->GetParent();
//    
////    ((TGWindow*)(f1->GetParent()))->MoveResize(0          ,0            , 0.66*width,heigth);
//    ((TGWindow*)(f2->GetParent()))->MoveResize(0.66*width, 0            , 0.33*width,0.33*heigth);
//    ((TGWindow*)(f3->GetParent()))->MoveResize(0.66*width, 0.33*heigth , 0.33*width,0.33*heigth);
//    ((TGWindow*)(f4->GetParent()))->MoveResize(0.66*width, 0.66*heigth , 0.33*width,0.33*heigth);

    
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();
    gEve->Redraw3D(kTRUE);

//    man->SetESDcolors(colors);
//    man->SetESDwidths(widths);
    man->SetESDdashNoRefit(dashNoRefit);
    man->SetESDdrawNoRefit(drawNoRefit);
    
    man->SetSaveViews(saveViews);
    man->SetESDtracksByCategory(false);
    man->SetESDtracksByType(true);
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