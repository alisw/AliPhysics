//
//  AliEveInit.cpp
//  xAliRoot
//
//  Created by Jeremi Niedziela on 01/06/15.
//
//

#include <AliEveInit.h>
#include <AliEveTrackCounter.h>
#include <AliEveEventManagerEditor.h>
#include <AliEveMultiView.h>
#include <AliEveMacroExecutor.h>
#include <AliEveMacro.h>
#include <AliEveGeomGentle.h>
#include <AliEveDataSourceOffline.h>
#include <AliEveDataSourceHLTZMQ.h>
#include <AliEveEventManager.h>

#include <AliCDBManager.h>

#include <TGrid.h>
#include <TROOT.h>
#include <TInterpreter.h>
#include <TEveWindowManager.h>
#include <TGTab.h>
#include <TPRegexp.h>
#include <TFolder.h>
#include <TSystemDirectory.h>
#include <TGFileBrowser.h>
#include <TEveMacro.h>
#include <TEveBrowser.h>

#include <iostream>

using namespace std;

AliEveInit::AliEveInit(const TString& path ,AliEveEventManager::EDataSource defaultDataSource,bool storageManager) :
fPath(path)
{
    //==============================================================================
    // Reading preferences from config file
    //==============================================================================
    
    TEnv settings;
    GetConfig(&settings);
    
    bool autoloadEvents   = settings.GetValue("events.autoload.set",false);   // set autoload by default
    bool fullscreen       = settings.GetValue("fullscreen.mode",false);       // hide left and bottom tabs
    
    TString ocdbStorage   = settings.GetValue("OCDB.default.path","local://$ALICE_ROOT/../src/OCDB");// default path to OCDB
    
    const int nDetectors = 2;
    const char* detectors[nDetectors] = {"ACO","MCH"};
    
    Info("AliEveInit",Form("\n\nOCDB path:%s\n\n",ocdbStorage.Data()));
    
    //==============================================================================
    // Event Manager and different data sources
    //==============================================================================
    
    AliEveEventManager *man = new AliEveEventManager(defaultDataSource);
    
    AliEveEventManager::SetCdbUri(ocdbStorage);
    
    if (gSystem->Getenv("ALICE_ROOT") != 0)
    {
        gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
        gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
    }
    
    AliEveDataSourceOffline *dataSourceOffline  = (AliEveDataSourceOffline*)man->GetDataSourceOffline();
    AliEveDataSourceHLTZMQ  *dataSourceHLT      = (AliEveDataSourceHLTZMQ*) man->GetDataSourceHLTZMQ();
    
    ImportMacros();
    Init();
    
    TEveUtil::AssertMacro("VizDB_scan.C");
    
    TEveBrowser *browser = gEve->GetBrowser();
    browser->ShowCloseTab(kFALSE);
    
    //==============================================================================
    // Geometry, scenes, projections and viewers
    //==============================================================================
    
    AliEveMultiView *mv = new AliEveMultiView();
    AliEveGeomGentle *geomGentle = new AliEveGeomGentle();
    
    mv->InitSimpleGeom(geomGentle->GetGeomGentleRphi(),true,false); // special geometry for RPhi projection
    mv->InitSimpleGeom(geomGentle->GetGeomGentle(),false,true);     // to be replaced by per-detector geometries
    mv->InitSimpleGeom(geomGentle->GetSimpleGeom("EMC"));
    mv->InitSimpleGeom(geomGentle->GetSimpleGeom("ACO"));
    mv->InitSimpleGeom(geomGentle->GetSimpleGeom("TRD"));
    
    if(settings.GetValue("MUON.show", true)){
        mv->InitSimpleGeom(geomGentle->GetSimpleGeom("MCH"),false);
    }
    
    AddMacros();
    
    //==============================================================================
    // Additional GUI components
    //==============================================================================
    
    TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    
    // QA viewer
    /*
     slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
     slot->StartEmbedding();
     new AliQAHistViewer(gClient->GetRoot(), 600, 400, kTRUE);
     slot->StopEmbedding("QA histograms");
     */
    //    browser->GetTabRight()->SetTab(1);
    browser->StartEmbedding(TRootBrowser::kBottom);
    new AliEveEventManagerWindow(man,storageManager,defaultDataSource);
    browser->StopEmbedding("EventCtrl");
    
    slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    TEveWindowTab *store_tab = slot->MakeTab();
    store_tab->SetElementNameTitle("WindowStore",
                                   "Undocked windows whose previous container is not known\n"
                                   "are placed here when the main-frame is closed.");
    gEve->GetWindowManager()->SetDefaultContainer(store_tab);
    
    
    //==============================================================================
    // AliEve objects - global tools
    //==============================================================================
    
    AliEveTrackCounter* g_trkcnt = new AliEveTrackCounter("Primary Counter");
    gEve->AddToListTree(g_trkcnt, kFALSE);
    
    
    //==============================================================================
    // Final stuff
    //==============================================================================
    
    
    // A refresh to show proper window.
    //    gEve->GetViewers()->SwitchColorSet();
    
    browser->MoveResize(0, 0, gClient->GetDisplayWidth(),gClient->GetDisplayHeight() - 32);
    gEve->Redraw3D(true);
    gSystem->ProcessEvents();
    
//    man->GotoEvent(0);
    
    gEve->EditElement(g_trkcnt);
    gEve->Redraw3D();
    
    // move and rotate sub-views
    browser->GetTabRight()->SetTab(1);
    TGLViewer *glv1 = mv->Get3DView()->GetGLViewer();
    TGLViewer *glv2 = mv->GetRPhiView()->GetGLViewer();
    TGLViewer *glv3 = mv->GetRhoZView()->GetGLViewer();
    
    glv1->CurrentCamera().RotateRad(-0.4, 1.0);
    glv2->CurrentCamera().Dolly(1, kFALSE, kFALSE);
    glv3->CurrentCamera().Dolly(1, kFALSE, kFALSE);
    
    // Fullscreen
    if(fullscreen){
        ((TGWindow*)gEve->GetBrowser()->GetTabLeft()->GetParent())->Resize(1,0);
        ((TGWindow*)gEve->GetBrowser()->GetTabBottom()->GetParent())->Resize(0,1);
        gEve->GetBrowser()->Layout();
    }
    
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();
    gEve->Redraw3D(true);
    
    man->SetAutoLoad(autoloadEvents);// set autoload by default
    
    if(defaultDataSource == AliEveEventManager::kSourceOffline){
        if(settings.GetValue("momentum.histograms.all.events.show",false))
        {
            ((AliEveDataSourceOffline*)man->GetDataSourceOffline())->GotoEvent(0);
            man->GetMomentumHistogramsDrawer()->DrawAllEvents();
        }
    }
}

void AliEveInit::Init()
{
    Info("AliEveInit","Adding standard macros");
    
    AliEveDataSourceOffline *dataSource = (AliEveDataSourceOffline*)AliEveEventManager::GetMaster()->GetDataSourceOffline();
    
    // Open event
    if (fPath.BeginsWith("alien:"))
    {
        if (gGrid != 0)
        {
            Info("AliEveInit::Init()", "TGrid already initializied. Skiping checks and initialization.");
        }
        else
        {
            Info("AliEveInit::Init()", "AliEn requested - connecting.");
            if (gSystem->Getenv("GSHELL_ROOT") == 0)
            {
                Error("AliEveInit::Init()", "AliEn environment not initialized. Aborting.");
                gSystem->Exit(1);
            }
            if (TGrid::Connect("alien") == 0)
            {
                Error("AliEveInit::Init()", "TGrid::Connect() failed. Aborting.");
                gSystem->Exit(1);
            }
        }
    }
    cout<<"Opening event -1 from "<<fPath.Data()<<endl;
    gEve->AddEvent(AliEveEventManager::GetMaster());
}

void AliEveInit::AddMacros()
{
    //==============================================================================
    // Registration of per-event macros
    //==============================================================================
    
    TEnv settings;
    GetConfig(&settings);
    
    bool showMuon         = settings.GetValue("MUON.show", true);              // show MUON's geom
    bool showEMCal        = settings.GetValue("EMCal.show", false);            // show EMCal and PHOS histograms
    bool drawClusters     = settings.GetValue("clusters.show",false);          // show clusters
    bool drawRawData      = settings.GetValue("rawData.show",false);           // show raw data
    bool drawHits         = settings.GetValue("hits.show",false);              // show hits
    bool drawDigits       = settings.GetValue("digits.show",false);            // show digits
    bool drawAD           = settings.GetValue("AD.show",false);                // show AD hits
    
    AliEveMacroExecutor *exec = AliEveEventManager::GetMaster()->GetExecutor();
    exec->RemoveMacros(); // remove all old macros
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Track",   "kine_tracks.C", "kine_tracks", "", kFALSE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits ITS", "its_hits.C",    "its_hits",    "", drawHits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits TPC", "tpc_hits.C",    "tpc_hits",    "", drawHits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits T0",  "t0_hits.C",     "t0_hits",     "", drawHits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits FMD", "fmd_hits.C",    "fmd_hits",    "", drawHits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits ACORDE", "acorde_hits.C",    "acorde_hits",    "", drawHits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits EMCAL", "emcal_hits.C",    "emcal_hits",    "", drawHits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits TOF",  "tof_hits.C",     "tof_hits",     "", drawHits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits TRD", "trd_hits.C",    "trd_hits",    "", drawHits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits VZERO", "vzero_hits.C",    "vzero_hits",    "", drawHits));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG ITS",     "its_digits.C",  "its_digits",  "", drawDigits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG TPC",     "tpc_digits.C",  "tpc_digits",  "", drawDigits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG TOF",     "tof_digits.C",  "tof_digits",  "", drawDigits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG HMPID",   "hmpid_digits.C","hmpid_digits","", drawDigits));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG FMD",     "fmd_digits.C",  "fmd_digits",  "", drawDigits));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ITS",     "its_raw.C",     "its_raw",     "", drawRawData));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TPC",     "tpc_raw.C",     "tpc_raw",     "", drawRawData));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TOF",     "tof_raw.C",     "tof_raw",     "", drawRawData));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW HMPID",   "hmpid_raw.C",   "hmpid_raw",   "", drawRawData));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW T0",      "t0_raw.C",      "t0_raw",      "", drawRawData));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW FMD",     "fmd_raw.C",     "fmd_raw",     "", drawRawData));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW VZERO",   "vzero_raw.C",   "vzero_raw",   "", drawRawData));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ACORDE",  "acorde_raw.C",  "acorde_raw",  "", drawRawData));

    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC ZDC",      "esd_zdc.C", "esd_zdc", "", kFALSE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters ITS", "its_clusters.C", "its_clusters","",     drawClusters));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters TPC", "tpc_clusters.C", "tpc_clusters","",     drawClusters));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters TRD", "trd_clusters.C", "trd_clusters","",     drawClusters));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters TOF", "tof_clusters.C", "tof_clusters","",     drawClusters));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters HMPID", "hmpid_clusters.C","hmpid_clusters","",drawClusters));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters PHOS", "phos_clusters.C","phos_clusters","",   drawClusters));
    
    if (showMuon)
    {
        exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM TrackRef MUON", "muon_trackRefs.C", "muon_trackRefs", "kTRUE", kFALSE));
        exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW MUON", "muon_raw.C", "muon_raw", "", drawRawData));
        exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG MUON", "muon_digits.C", "muon_digits", "", drawDigits));
        exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters MUON", "muon_clusters.C", "muon_clusters", "", drawClusters));
    }
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD AD", "ad_esd.C", "ad_esd", "", drawAD));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD EMCal", "emcal_esdclustercells.C", "emcal_esdclustercells", "", showEMCal));
}

void AliEveInit::ImportMacros()
{
    // Put macros in the list of browsables, add a macro browser to
    // top-level GUI.
    
    TString  hack = gSystem->pwd(); // Problem with TGFileBrowser cding
    
    TString macdir("$(ALICE_ROOT)/EVE/macros/data");
    gSystem->ExpandPathName(macdir);
    macdir = "$(ALICE_ROOT)/EVE/macros/common";
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
    gSystem->cd(hack);
}

void AliEveInit::GetConfig(TEnv *settings)
{
    TEveException kEH("AliEveInit::GetConfig");
    
    if(settings->ReadFile(Form("%s/eve_config",gSystem->Getenv("HOME")), kEnvUser) < 0)
    {
        Warning(kEH," could not find eve_config in home directory! Trying in $ALICE_ROOT/EVE/EveBase/");
        if(settings->ReadFile(Form("%s/EVE/EveBase/eve_config",gSystem->Getenv("ALICE_ROOT")), kEnvUser) < 0)
        {
            Error(kEH,"could not find eve_config file!.");
            exit(0);
        }
        else{
            Info(kEH,"Read config from standard location");
        }
    }
    else{
        Info(kEH,"Read config from home directory");
    }
}



